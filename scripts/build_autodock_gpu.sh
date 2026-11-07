#!/bin/bash
# Build AutoDock-GPU (CUDA) with a higher limit on ligand atom types.
#
# The stock build supports 32 atom types (MAX_NUM_OF_ATYPES in common/defines.h). These
# are shared by the anchor atoms and the ordinary types of the ligand's other atoms, so a ligand with a very large
# shared core can exceed them. This builds a copy of the source with a higher limit; the source folder you point
# it at is not changed.
#
# Usage:  scripts/build_autodock_gpu.sh <AutoDock-GPU source folder> [max atom types] [output folder]
#   max atom types   default 64
#   output folder    default <repo>/bin/adgpu_atypes<N>
#
# Environment variables:
#   TARGETS   CUDA compute capabilities to build for, in one binary (default "80 90": the A100 GPUs of the
#             ga100 partition and the H100 GPUs of gh100). A binary only runs on the capabilities it contains.
#   NUMWI     work-items per block: 32, 64, 128 or 256 (default 128)
#   JOBS      parallel compile jobs (default 4)
# CUDA_MODULES, GCC_MODULE   modules to load (defaults: CUDA/12.1.1 and GCC/12.3.0; CUDA 12.1 does not accept GCC 13)
#
# The output folder gets the binary (autodock_gpu_<NUMWI>wi), a file max_atom_types with the limit and a file
# targets with the compute capabilities, which dock_cli.py reads. Point the docking at it with --adg_path.

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOURCE="${1:?Usage: $0 <AutoDock-GPU source folder> [max atom types] [output folder]}"
MAX_TYPES="${2:-64}"
OUTPUT="${3:-$REPO_DIR/bin/adgpu_atypes$MAX_TYPES}"
TARGETS="${TARGETS:-80 90}"
NUMWI="${NUMWI:-128}"
JOBS="${JOBS:-4}"
CUDA_MODULES="${CUDA_MODULES:-CUDA/12.1.1}"
GCC_MODULE="${GCC_MODULE:-GCC/12.3.0}"

if ! [[ "$MAX_TYPES" =~ ^[0-9]+$ ]] || [ "$MAX_TYPES" -lt 32 ] || [ "$MAX_TYPES" -gt 256 ]; then
    echo "Max atom types must be a whole number from 32 to 256 (got '$MAX_TYPES')." >&2
    exit 1
fi
if [ ! -f "$SOURCE/common/defines.h" ] || [ ! -f "$SOURCE/Makefile" ]; then
    echo "$SOURCE does not look like an AutoDock-GPU source folder (common/defines.h or Makefile is missing)." >&2
    exit 1
fi

if type module >/dev/null 2>&1; then
    module load "$GCC_MODULE" $CUDA_MODULES
fi
if ! command -v nvcc >/dev/null 2>&1; then
    echo "nvcc was not found. Load a CUDA module first (for example: module load CUDA/12.1.1)." >&2
    exit 1
fi
CUDA_ROOT="${EBROOTCUDA:-$(dirname "$(dirname "$(command -v nvcc)")")}"
# CUDA 12 installs may keep the headers under targets/; libcuda is only a link-time stub here, the real
# driver library is found at run time on the GPU node.
CUDA_INCLUDE="$CUDA_ROOT/include"
[ -f "$CUDA_INCLUDE/cuda.h" ] || CUDA_INCLUDE="$CUDA_ROOT/targets/x86_64-linux/include"
[ -f "$CUDA_INCLUDE/cuda.h" ] || { echo "cuda.h was not found under $CUDA_ROOT." >&2; exit 1; }
CUDA_LIB="$CUDA_ROOT/lib64"
ls "$CUDA_LIB"/libcudart.so* >/dev/null 2>&1 || CUDA_LIB="$CUDA_ROOT/targets/x86_64-linux/lib"
ls "$CUDA_LIB"/libcudart.so* >/dev/null 2>&1 || { echo "libcudart was not found under $CUDA_ROOT." >&2; exit 1; }
[ -d "$CUDA_ROOT/stubs/lib64" ] && export LIBRARY_PATH="$CUDA_ROOT/stubs/lib64${LIBRARY_PATH:+:$LIBRARY_PATH}"

BUILD="$(mktemp -d "${TMPDIR:-/tmp}/adgpu_build.XXXXXX")"
trap 'rm -rf "$BUILD"' EXIT
cp -r "$SOURCE/." "$BUILD/"
rm -rf "$BUILD/.git" "$BUILD/bin"
mkdir -p "$BUILD/bin"

# Exactly one definition must change, otherwise this source differs from what the script expects.
if [ "$(grep -c '^#define MAX_NUM_OF_ATYPES ' "$BUILD/common/defines.h")" != "1" ]; then
    echo "Could not find a single 'MAX_NUM_OF_ATYPES' definition in common/defines.h." >&2
    exit 1
fi
sed -i "s/^#define MAX_NUM_OF_ATYPES .*/#define MAX_NUM_OF_ATYPES          $MAX_TYPES/" "$BUILD/common/defines.h"
grep -n '^#define MAX_NUM_OF_ATYPES ' "$BUILD/common/defines.h"

echo "Building for compute capability $TARGETS with $NUMWI work-items (CUDA in $CUDA_ROOT)..."
( cd "$BUILD" && make -j"$JOBS" DEVICE=CUDA NUMWI="$NUMWI" TARGETS="$TARGETS" \
    GPU_INCLUDE_PATH="$CUDA_INCLUDE" GPU_LIBRARY_PATH="$CUDA_LIB" )

BINARY="$BUILD/bin/autodock_gpu_${NUMWI}wi"
if [ ! -x "$BINARY" ]; then
    echo "The build did not produce $BINARY." >&2
    exit 1
fi
mkdir -p "$OUTPUT"
cp "$BINARY" "$OUTPUT/"
[ -x "$BUILD/bin/adgpu_analysis" ] && cp "$BUILD/bin/adgpu_analysis" "$OUTPUT/"   # the result analysis tool
echo "$MAX_TYPES" > "$OUTPUT/max_atom_types"
echo "$TARGETS" > "$OUTPUT/targets"
echo "Built $OUTPUT/autodock_gpu_${NUMWI}wi with a limit of $MAX_TYPES atom types for compute capability $TARGETS."
