# Anchored Docking Workflow

An automated pipeline for anchored docking using AutoDock-GPU. By overlaying ligands that share a common core with a template ligand, this serves as an alternative to core-constrained docking.

Based on the [AutoDock-GPU Anchored Docking Guide](https://github.com/ccsb-scripps/AutoDock-GPU/wiki/Anchored-docking).

> **Note:** This project is under active development and has not been experimentally validated.

New to the command line? Follow the step-by-step [Beginner's Guide](docs/BEGINNERS_GUIDE.md).

## Prerequisites

1. [Pixi](https://pixi.sh/) (0.83.0 or newer) — installs RDKit, Meeko, AutoGrid and ADFRsuite
2. [AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU/releases/) — the docking engine. On the HPC it is provided as a module and must be used together with GCC 13.2:

   ```bash
   module load GCC/13.2.0 AutoDock-GPU/1.5.3-CUDA
   ```

   Docking needs a GPU, e.g. `--partition=ga100 --gres=gpu:1`.

## Installation

```bash
git clone https://github.com/Chemical-Biology-STP/anchored_docking_workflow.git
cd anchored_docking_workflow
pixi install --all
```

This creates two local `.pixi` environments:

| Environment | Contents |
|---|---|
| `default` | Python 3.12, RDKit 2026.09, Meeko 0.8, AutoGrid 4.2.9 |
| `adfr` | ADFRsuite 1.0 (`prepare_receptor`) |

ADFRsuite is in its own environment because its conda package bundles Python 2.7 and would otherwise overwrite `python` in the default environment.

### Tool locations

Nothing needs to be configured by hand:

- `autogrid4` is taken from the pixi `default` environment.
- `prepare_receptor` is taken from `.pixi/envs/adfr` (override with `--adfr_path`).
- The AutoDock-GPU binary (`adgpu` or `autodock_gpu_*wi`) is found on `PATH` after loading the module (override with `--adg_path <dir>`).

### Scripts directory

The `scripts/` directory contains helper scripts required by the pipeline:

- `write-gpf.py` — from [diogomart/write-autogrid-config](https://github.com/diogomart/write-autogrid-config)
- `addbias.py` and `insert_type_in_fld.py` — see [AutoDock-GPU#283](https://github.com/ccsb-scripps/AutoDock-GPU/issues/283)

These are already included in this repository.

## Usage

Run on a GPU node (interactively or via `sbatch`, see below):

```bash
module load GCC/13.2.0 AutoDock-GPU/1.5.3-CUDA

pixi run python anchored_docking.py \
    --template_file test/crystal_ligand.sdf \
    --input_file test/CAT-13f.sdf \
    --receptor_file test/Bace1_protein.pdb \
    --docking_dir test/output
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
| `--adg_path` | search `PATH` | Directory containing the AutoDock-GPU binary (`adgpu` or `autodock_gpu_*wi`) |
| `--adfr_path` | `.pixi/envs/adfr` | ADFRsuite installation directory |
| `--scripts_path` | `scripts` | Path to the helper scripts directory |
| `--force` | `False` | Overwrite the output directory if it already exists |
| `--resume` | `False` | Continue a previous run in `--docking_dir`: ligands that already have a docked pose are skipped and the others are redone. Refused if the inputs or docking settings differ from the previous run. Cannot be combined with `--force` |
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
9. Checks the pose against the anchoring (see below) and generates a PyMOL script (`view_results.pml`) for visualising the results

### Pose check

Docking only pulls the anchor atoms towards their targets with a soft potential, so a pose can come out with the core in the wrong place and still be written as a normal result. After each ligand is docked, `pose_check.py` measures the top pose:

- **Anchor deviation:** the largest distance between a docked anchor atom and its target in the template.
- **Core RMSD:** the RMSD over all atoms of the shared core, against the template.

A pose is flagged (`CHECK POSE`) if an anchor is more than 2.0 Å from its target or the core RMSD is above 2.5 Å. Two real AutoDock-GPU runs of the example ligand gave anchor deviations of 1.2 to 1.3 Å and core RMSDs of 0.8 Å, which is why the anchor limit is 2.0 Å; both limits are judgements, not calibrated values, so look at the numbers as well as the flag. The thresholds are constants at the top of `pose_check.py`. For a symmetric core the mapping that fits best is used, so the numbers are the most favourable reading of the pose.

### Failures and resuming

A problem with one ligand (it cannot be read, has no shared core, fails 3D embedding, or an external program fails) is recorded and the remaining ligands are still docked. The script then exits with a non-zero code so the job shows as failed. Per-ligand status is written to `<docking_dir>/summary.tsv`, and the reason for a failure to `<docking_dir>/<index>/error.txt`. Fix the cause and re-run with `--resume` to redo only the ligands that have no docked pose.

A ligand whose 3D coordinates cannot be generated is reported as failed. Earlier versions carried on with a flat 2D geometry, which Meeko accepted without complaint. Each embedding attempt is limited to 30 seconds.

## Viewing results

After docking completes, each output subdirectory contains a `view_results.pml` script. Open it in [PyMOL](https://pymol.org/) to see the receptor, template ligand, and docked pose together:

```bash
pymol test/output/0/view_results.pml
```

This loads all three structures with sensible colours and styling, and saves a `docking_results.pse` session file you can reopen later without the script.

## Interactive CLI

The easiest way to run the workflow is the interactive command-line interface:

```bash
pixi run dock
```

(or `pixi shell` followed by `dock`). It asks for the receptor, template ligand and ligand files one at a time, with Tab completion for paths, and checks each input as you go:

- **Receptor:** warns if it has no hydrogens, still contains other molecules, or has several models or alternate locations.
- **Template:** checks it has real 3D coordinates and sits in the receptor, which catches a template and receptor in different coordinate systems. Warns if the file holds several molecules (only the first is used).
- **Ligands:** runs the same core and anchor search the docking uses on every ligand, tries the same 3D embedding, and flags ligands with no or a tiny shared core, no 3D geometry, several fragments (salts), or atoms AutoDock has no atom type for. Flagged ligands can be left out before the job starts.
- **Earlier runs:** if the output folder holds a previous run with identical inputs and settings, you can resume it so finished ligands are skipped. If anything differs, resuming is refused and the reason is shown.
- **Results:** the final table shows each ligand's status, the pose check numbers, and the reason for any failure, instead of only whether a file exists.
- **Settings:** the anchor count (at most 17, from AutoDock-GPU's 32 atom types), box size (at most 95 Å, from its 256 grid points), random seed (0 or higher) and time limit are range-checked. Manual core patterns are validated against the template.
- **Resources:** warns when the estimated disk use is large, refuses unwritable or dangerous output folders, and requires the ADFRsuite environment (`pixi install --all`) before starting.

It then asks how to run:

1. **Submit to SLURM:** writes a job script (GPU partition `ga100`, modules `GCC/13.2.0` and `AutoDock-GPU/1.5.3-CUDA`, a time limit scaled to the number of ligands), submits it with `sbatch`, and can watch the job and summarise the results when it finishes.
2. **Run now:** runs the docking in the current terminal. This needs a GPU and the AutoDock-GPU module already loaded.
3. **Save only:** writes the job script without running anything.

Nothing is submitted until you confirm on the final review screen. The output folder is emptied if it already exists, so the CLI refuses a folder that contains your input files or this repository and asks before deleting a non-empty one.

The CLI is a front end for `anchored_docking.py`; you can still call that script directly with the flags listed above.

### Tests

The CLI's checks have automated tests that need no GPU or SLURM:

```bash
pixi run -e dev test
```

The limits above come from the compiled AutoDock-GPU 1.5.3 build used on the cluster (`common/defines.h`). If you switch to a different build, update the constants at the top of `dock_cli.py`.

## Example

Using the included Bace1 test case:

```bash
module load GCC/13.2.0 AutoDock-GPU/1.5.3-CUDA

pixi run python anchored_docking.py \
    --template_file test/crystal_ligand.sdf \
    --input_file test/CAT-13f.sdf \
    --receptor_file test/Bace1_protein.pdb \
    --docking_dir test/output \
    --force
```

<img src="docked_pose.png" width="60%">

- Green: receptor (Bace1)
- Cyan: template ligand (CAT-4j)
- Magenta: docked CAT-13f
- Yellow: CAT-13f pose from the [Wang-FEP-dataset](https://github.com/ohuelab/Wang-FEP-dataset)

## Acknowledgements

Refactored and improved by Daniella Hares and Yew Mun Yip from the Chemical Biology STP.
