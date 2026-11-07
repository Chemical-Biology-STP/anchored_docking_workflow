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
| `--anchor_num` | `3` | Number of anchor points to use. In `spread` mode it is a maximum: a ligand whose core has fewer usable atoms uses all of them. `0` uses every usable atom, which fails for a ligand that would need more atom types than AutoDock-GPU supports |
| `--anchor_mode` | `random` | Anchor selection mode: `random`, `area` (always 3 anchors, maximising triangle area) or `spread` (up to `--anchor_num` core atoms that are as far apart from each other as possible; deterministic, no seed needed). The `pixi run dock` wizard recommends `spread` with the atom-type limit as the maximum, so every usable core atom is used where AutoDock-GPU allows it |
| `--query_smarts` | `None` | Manually specify a SMARTS pattern for the common substructure instead of using MCS |
| `--anchor_indices` | `None` | Specify exact anchor atom indices (0-based) within the query SMARTS. Requires `--query_smarts` |
| `--size` | `20.0` | Size of the docking box in Angstroms |
| `--samples` | `100` | Number of random samples when using `area` anchor mode |
| `--random_seed` | `None` | Random seed for reproducible anchor selection |
| `--adg_path` | search `PATH` | Directory containing the AutoDock-GPU binary (`adgpu` or `autodock_gpu_*wi`) |
| `--adfr_path` | `.pixi/envs/adfr` | ADFRsuite installation directory |
| `--scripts_path` | `scripts` | Path to the helper scripts directory |
| `--force` | `False` | Overwrite the output directory if it already exists |
| `--template_view_file` | `None` | Display only: the template with its hydrogens, for a heavy-atom template whose bond orders are unknown (the CLI sets this for templates taken from a structure file). Without it, missing template hydrogens are added with RDKit. Never used for docking |
| `--add_hydrogens` | `False` | Let ADFR add hydrogens to the receptor if it has none (`prepare_receptor -A checkhydrogens`). Receptors that already have hydrogens are unchanged. Counts as a docking setting for `--resume` |
| `--resume` | `False` | Continue a previous run in `--docking_dir`: ligands that already have a docked pose are skipped and the others are redone. Refused if the inputs or docking settings differ from the previous run. Cannot be combined with `--force` |
| `--clean` | `False` | Remove intermediate files after docking |
| `--verbose` | `False` | Print additional debug information (SMILES, SMARTS, etc.) |

## How it works

1. Finds the maximum common substructure (MCS) between the template and input ligand (or uses a user-provided SMARTS)
2. Selects anchor atoms — unique, non-equivalent heavy atoms in the common core (atoms with a symmetric twin are never used, because their mapping to the template is ambiguous). Each anchor needs its own atom type, and AutoDock-GPU 1.5.3 supports 32 atom types per ligand in total, counting the anchors **and** the ordinary types of the ligand's other atoms. How many anchors fit therefore depends on the ligand (22 for a 26-atom example, 29 for a 46-atom one). In `spread` mode the script lowers the number of anchors until the ligand fits and says so in the results; in the other modes such a ligand fails with an explanation. A build with a higher limit lifts this (see [below](#autodock-gpu-with-a-higher-atom-type-limit))
3. Generates custom atom types for anchor atoms and writes a Meeko parameters file
4. Prepares the ligand (PDBQT via Meeko) and receptor (PDBQT via ADFR)
5. Runs AutoGrid to generate grid maps
6. Adds Gaussian bias potentials to grid maps at anchor positions
7. Runs AutoDock-GPU with the biased maps
8. Exports the docked pose as SDF
9. Checks the pose against the anchoring (see below), then writes the PyMOL scripts, the combined result tables and the versions of the protein and template with every hydrogen (see [Output](#output))

### Pose check

Docking only pulls the anchor atoms towards their targets with a soft potential, so a pose can come out with the core in the wrong place and still be written as a normal result. After each ligand is docked, `pose_check.py` measures the top pose:

- **Anchor deviation:** the largest distance between a docked anchor atom and its target in the template.
- **Core RMSD:** the RMSD over all atoms of the shared core, against the template.

A pose is flagged (`CHECK POSE`) if an anchor is more than 2.0 Å from its target or the core RMSD is above 2.5 Å. Two real AutoDock-GPU runs of the example ligand gave anchor deviations of 1.2 to 1.3 Å and core RMSDs of 0.8 Å, which is why the anchor limit is 2.0 Å; both limits are judgements, not calibrated values, so look at the numbers as well as the flag. The thresholds are constants at the top of `pose_check.py`. For a symmetric core the mapping that fits best is used, so the numbers are the most favourable reading of the pose.

### Failures and resuming

A problem with one ligand (it cannot be read, has no shared core, fails 3D embedding, or an external program fails) is recorded and the remaining ligands are still docked. The script then exits with a non-zero code so the job shows as failed. Per-ligand status is written to `<docking_dir>/summary.tsv`, and the reason for a failure to `<docking_dir>/<index>/error.txt`. Fix the cause and re-run with `--resume` to redo only the ligands that have no docked pose.

A ligand whose 3D coordinates cannot be generated is reported as failed. Earlier versions carried on with a flat 2D geometry, which Meeko accepted without complaint. Each embedding attempt is limited to 30 seconds.

## Protein preparation

`pixi run dock` prepares the protein before docking (`protein_prep.py`); `anchored_docking.py` itself only uses the receptor it is given.

- **Gaps and missing atoms** are reported for every receptor: consecutive residues whose C–N distance exceeds 2.0 Å (a real peptide bond is 1.33 Å), and standard residues lacking heavy atoms. They are shown as warnings when within 10 Å (gap) or 6 Å (residue) of the template.
- **Protonation** is offered when the receptor has no hydrogens, at a pH you choose (default **7.4**), using [PDB2PQR](https://pdb2pqr.readthedocs.io) with [PROPKA](https://github.com/jensengroup/propka): predicted pKa values set the form of each histidine (HID/HIE/HIP), aspartate, glutamate, cysteine (disulfides become CYX), lysine and tyrosine, hydrogens are placed by hydrogen-bond-network optimisation, and missing heavy atoms are rebuilt. Receptors that already have hydrogens are used unchanged.
- **Heavy atoms can move.** PDB2PQR flips Asn, Gln and His side chains when that improves hydrogen bonds (on the 1M17 test protein, two residues changed; the largest shift in a plain Bace1 run was 2.4 Å). Flipped residues are listed in the report.
- The result is saved as `<job>_receptor.pdb` with `<job>_protonation_report.txt` (states, pKa values near the template, flips, rebuilt atoms, gaps). If PDB2PQR fails, or the user declines, the CLI offers ADFR's rule-based hydrogens (`--add_hydrogens`), which make every histidine cationic.
- A kept ion or cofactor passes through PDB2PQR unchanged and gets no hydrogens.

## AutoDock-GPU with a higher atom-type limit

The AutoDock-GPU module here is compiled for 32 atom types per ligand (`MAX_NUM_OF_ATYPES` in `common/defines.h`). A ligand needs one type per anchor plus the types of its other atoms, so a core of about 30 usable atoms can exceed it. The stock binary then stops with `Too many ligand atom types (more than 32)`. (AutoDock-GPU itself exits with code 0 after that error; the script now catches it.)

`scripts/build_autodock_gpu.sh` builds a copy of an AutoDock-GPU source tree with a higher limit and leaves the source untouched:

```bash
scripts/build_autodock_gpu.sh /path/to/AutoDock-GPU 64
```

- It needs the CUDA and GCC modules. The defaults are `CUDA/12.1.1` and `GCC/12.3.0`, because CUDA 12.1 rejects GCC 13. It takes a few minutes and does not need a GPU.
- It builds one binary for compute capabilities 80 and 90 (the A100 GPUs of `ga100` and the H100 GPUs of `gh100`) unless `TARGETS` is set, with 128 work-items like the module. A binary only runs on the capabilities it contains: an 80-only build fails on an H100 with `Assertion '0' failed` in `SetKernelsGpuData`. The result is in `bin/adgpu_atypes64/`, together with `max_atom_types` and `targets` files; the wizard reads the latter and refuses a build that does not match the partition's GPUs.
- The cluster module's own binary contains code for compute capabilities 50 to 90, so it needs no such care.
- **As a cluster module:** `AutoDock-GPU/1.5.3-CUDA-atypes64` is installed next to the original module, built by this script from the same source tree. It loads the same GCC and CUDA modules, contains code for A100 and H100 GPUs, and also provides `adgpu_analysis`. The original `AutoDock-GPU/1.5.3-CUDA` is unchanged and stays the default for an unversioned `module load AutoDock-GPU` (a `.modulerc` records that). `pixi run dock` finds the module's binaries in the shared software folder as well as in `bin/`.
- **The wizard chooses it by itself.** For every ligand whose core is large enough that it could matter (more than 17 usable anchor atoms, since the other atoms of a ligand have at most 15 different types and 17 + 15 = 32), `pixi run dock` runs the same anchor choice and Meeko preparation as the job and counts the atom types exactly. If any ligand needs more than the 32 of the cluster's default AutoDock-GPU, it uses the higher-limit build without asking, and says which ligands and why. The review screen shows which AutoDock-GPU the job uses, and the advanced settings let you go back to the default. If no such build exists, it warns that those ligands will get fewer anchors; if a ligand cannot be measured, it says so.
- Use a local build by hand with `--adg_path bin/adgpu_atypes64`. The script reads the limit from the `max_atom_types` file beside the binary.
- Tested on an A100 node and on an H100 node (`gh100`, reservation `h100`): a 30-anchor ligand gets 29 anchors with the module's binary and all 30 with the 64-type build, on both GPU types. Nothing else was changed in the AutoDock-GPU source, and for a ligand that fits 32 types the two binaries gave the same result in the one comparison made (free energy −9.32 against −9.33 kcal/mol, same anchors, same pose check).

## Output

The job writes no PyMOL session itself: it writes PyMOL **scripts**, and PyMOL saves the session when you run one.

In the output folder:

| File | What it is |
|---|---|
| `results.csv` | One row per ligand: status, free energy and its inter/intra parts, cluster size, number of poses, anchor deviation, core RMSD, failure reason. |
| `top_poses.sdf` | The best pose of every docked ligand, with the scores and pose check as data fields. |
| `all_poses.sdf` | Every cluster-lead pose of every docked ligand (`pose_rank` 1, 2, ...). |
| `view_all.pml` | Opens the protein, the template and all docked ligands in one PyMOL session (flagged poses in orange) and saves `docking_all.pse`. |
| `receptor_with_hydrogens.pdb`, `template_with_hydrogens.sdf` | The protein and template used for viewing, with every hydrogen. |
| `summary.tsv`, `run_settings.json` | Per-ligand status for the run, and the settings used (for `--resume`). |
| `<index>/ligand_docked.sdf` | That ligand's cluster-lead poses, complete molecules with scores in a `meeko` property. |
| `<index>/view_results.pml` | Opens just that ligand with the protein and template, and saves `docking_results.pse`. |
| `<index>/pose_check.json`, `anchor_positions.json`, `error.txt` | Pose check numbers, anchor targets, and the failure details if the ligand failed. |
| `<index>/ligand.dlg` etc. | The raw AutoDock-GPU output and intermediates (removed with `--clean`). |

Run the PyMOL scripts from inside the folder that holds them:

```bash
cd test/output && pymol view_all.pml        # everything in one session
cd test/output/0 && pymol view_results.pml  # one ligand
```

The combined files can be rebuilt from an existing run with `pixi run python results.py <docking_dir>`.

### Hydrogens

AutoDock merges non-polar hydrogens into their carbons, so its own files are not what you want to look at. For viewing, every molecule has all its hydrogens:

- **Ligands:** Meeko already writes the complete molecule, so the poses in `ligand_docked.sdf` have every hydrogen (checked: 23 of 23 for the example ligand, bond lengths 1.04 to 1.10 Å).
- **Protein:** if your receptor file had hydrogens it is shown as supplied. If ADFR added them (`--add_hydrogens`), ADFR is run a second time without merging the non-polar ones. The polar hydrogens come out at exactly the coordinates used for docking (checked: 538 of 538 at 0.0000 Å), plus the non-polar ones.
- **Template:** shown with its hydrogens. A template given as SDF gets any missing hydrogens added by RDKit. A template taken from a structure file has no bond orders, so its hydrogens come from the PDB's chemical component dictionary for that residue name (fetched from files.rcsb.org when you run the CLI). This is checked: every heavy atom must match and the hydrogen count must equal the dictionary's, otherwise it is shown without hydrogens and the reason is printed. The docking itself always uses the heavy-atom template.

The protein's non-polar hydrogens and its histidine protonation are a quick rule-based placement, not an optimised one.

## Interactive CLI

The easiest way to run the workflow is the interactive command-line interface:

```bash
pixi run dock
```

(or `pixi shell` followed by `dock`). It starts from either a **protein-ligand complex** (PDB or mmCIF, for example straight from rcsb.org) or a cleaned receptor plus a separate template SDF. From a complex it does the clean-up for you (see below). It asks for each file in turn, with Tab completion for paths, and checks each input as you go:

- **Complex:** finds the molecules in the file and asks which is the template (suggesting the most likely one), asks which protein chains to keep (suggesting those near the template), and asks about each ion, cofactor or other molecule near the template. Waters, crystallisation additives (sulfate, glycerol, PEG, ...), sugars, other copies of the template and molecules more than 12 Å away are removed without asking; everything removed is listed. It keeps the first model and the first alternate location, writes a cleaned receptor and the extracted template next to the results (your input is never modified), and offers to add hydrogens if there are none (ADFR's `-A checkhydrogens`, a quick guess that does not optimise histidine protonation). The template's bond orders are not in a structure file, so it is written with single bonds; this does not affect the shared core or the anchors. Only the first model of mmCIF/PDB files is used and nucleic acid chains are dropped.
- **Receptor:** warns if it has no hydrogens, still contains other molecules, or has several models or alternate locations.
- **Template:** checks it has real 3D coordinates and sits in the receptor, which catches a template and receptor in different coordinate systems. Warns if the file holds several molecules (only the first is used).
- **Ligands:** runs the same core and anchor search the docking uses on every ligand, tries the same 3D embedding, and flags ligands with no or a tiny shared core, no 3D geometry, several fragments (salts), or atoms AutoDock has no atom type for. Flagged ligands can be left out before the job starts.
- **Earlier runs:** if the output folder holds a previous run with identical inputs and settings, you can resume it so finished ligands are skipped. If anything differs, resuming is refused and the reason is shown.
- **Results:** the final table shows each ligand's status, the pose check numbers, and the reason for any failure, instead of only whether a file exists.
- **Settings:** the anchor count (limited by AutoDock-GPU's 32 atom types, shared with the ligand's other atom types), box size (at most 95 Å, from its 256 grid points), random seed (0 or higher) and time limit are range-checked. Manual core patterns are validated against the template.
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
