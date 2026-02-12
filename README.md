# Anchored Docking Workflow

An automated pipeline for anchored docking using AutoDock-GPU. By overlaying ligands that share a common core with a template ligand, this serves as an alternative to core-constrained docking.

Based on the [AutoDock-GPU Anchored Docking Guide](https://github.com/ccsb-scripps/AutoDock-GPU/wiki/Anchored-docking).

> **Note:** This project is under active development and has not been experimentally validated.

## Prerequisites

You need the following installed and accessible on your system:

1. [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU/releases/) — the docking engine (`autodock_gpu_128wi` or similar binary)
2. [AutoGrid4](https://autodock.scripps.edu/download-autodock4/) — for grid map generation
3. [ADFRsuite](https://ccsb.scripps.edu/adfr/downloads/) — for receptor preparation (`prepare_receptor`)
4. [Pixi](https://pixi.sh/) — package manager for Python dependencies

## Installation

```bash
git clone https://github.com/Chemical-Biology-STP/anchored_docking_workflow.git
cd anchored_docking_workflow
pixi install
```

This installs the Python dependencies (RDKit, Meeko 0.6.1) into a local `.pixi` environment.

### Setting up binaries

The script expects paths to AutoDock-GPU and ADFRsuite. You can either:

- Pass them via `--adg_path` and `--adfr_path` flags, or
- Create symlinks in a `bin/` directory:

```bash
mkdir -p bin
ln -s /path/to/autodock_gpu_128wi bin/adgpu
ln -s /path/to/autogrid4 bin/autogrid4
```

### Scripts directory

The `scripts/` directory contains helper scripts required by the pipeline:

- `write-gpf.py` — from [diogomart/write-autogrid-config](https://github.com/diogomart/write-autogrid-config)
- `addbias.py` and `insert_type_in_fld.py` — see [AutoDock-GPU#283](https://github.com/ccsb-scripps/AutoDock-GPU/issues/283)

These are already included in this repository.

## Usage

```bash
pixi run dock
```

Or run directly:

```bash
pixi run python anchored_docking.py \
    --template_file test/crystal_ligand.sdf \
    --input_file test/CAT-13f.sdf \
    --receptor_file test/Bace1_protein.pdb \
    --docking_dir test/output \
    --adg_path /path/to/autodock-gpu/bin \
    --adfr_path /path/to/ADFRsuite
```

## Parameters

### Required

| Parameter | Description |
|---|---|
| `--template_file` | SDF file of the reference/crystal ligand used to define anchor positions |
| `--input_file` | SDF file of the ligand(s) to dock (can contain multiple molecules) |
| `--receptor_file` | PDB file of the receptor protein |
| `--docking_dir` | Output directory for docking results |

### Optional

| Parameter | Default | Description |
|---|---|---|
| `--anchor_num` | `3` | Number of anchor points to use |
| `--anchor_mode` | `random` | Anchor selection mode: `random` or `area` (area always uses 3 anchors, maximising triangle area) |
| `--query_smarts` | `None` | Manually specify a SMARTS pattern for the common substructure instead of using MCS |
| `--anchor_indices` | `None` | Specify exact anchor atom indices (0-based) within the query SMARTS. Requires `--query_smarts` |
| `--size` | `20.0` | Size of the docking box in Angstroms |
| `--samples` | `100` | Number of random samples when using `area` anchor mode |
| `--random_seed` | `None` | Random seed for reproducible anchor selection |
| `--adg_path` | `autodock-gpu` | Path to directory containing the AutoDock-GPU binary (`adgpu`) |
| `--adfr_path` | `ADFRsuite-1.0` | Path to ADFRsuite installation directory |
| `--scripts_path` | `scripts` | Path to the helper scripts directory |
| `--force` | `False` | Overwrite the output directory if it already exists |
| `--clean` | `False` | Remove intermediate files after docking |
| `--verbose` | `False` | Print additional debug information (SMILES, SMARTS, etc.) |

## How it works

1. Finds the maximum common substructure (MCS) between the template and input ligand (or uses a user-provided SMARTS)
2. Selects anchor atoms — unique, non-equivalent heavy atoms in the common core
3. Generates custom atom types for anchor atoms and writes a Meeko parameters file
4. Prepares the ligand (PDBQT via Meeko) and receptor (PDBQT via ADFR)
5. Runs AutoGrid to generate grid maps
6. Adds Gaussian bias potentials to grid maps at anchor positions
7. Runs AutoDock-GPU with the biased maps
8. Exports the docked pose as SDF
9. Generates a PyMOL script (`view_results.pml`) for visualising the results

## Viewing results

After docking completes, each output subdirectory contains a `view_results.pml` script. Open it in [PyMOL](https://pymol.org/) to see the receptor, template ligand, and docked pose together:

```bash
pymol test/output/0/view_results.pml
```

This loads all three structures with sensible colours and styling, and saves a `docking_results.pse` session file you can reopen later without the script.

## Generating a SLURM submission script

An interactive generator with tab completion is included:

```bash
pixi run python generate_slurm.py
```

It walks you through SLURM settings (partition, GPUs, time limit, modules) and all docking parameters, then writes a ready-to-submit script. Tab works for file paths and choice fields.

## Example

Using the included Bace1 test case:

```bash
pixi run python anchored_docking.py \
    --template_file test/crystal_ligand.sdf \
    --input_file test/CAT-13f.sdf \
    --receptor_file test/Bace1_protein.pdb \
    --docking_dir test/output \
    --force \
    --adg_path $PWD/bin \
    --adfr_path /path/to/ADFRsuite
```

<img src="docked_pose.png" width="60%">

- Green: receptor (Bace1)
- Cyan: template ligand (CAT-4j)
- Magenta: docked CAT-13f
- Yellow: CAT-13f pose from the [Wang-FEP-dataset](https://github.com/ohuelab/Wang-FEP-dataset)
