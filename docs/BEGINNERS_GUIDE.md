# Beginner's Guide: Anchored Docking on the HPC

This guide takes you from "I have a protein and some ligands" to "I have docked poses I can look at", assuming you have never used a command line before. Follow the steps in order and copy the commands exactly.

> **Status:** this workflow is under active development and has not been experimentally validated. Treat results as hypotheses to check, not as facts.

## Contents

1. [What this does](#1-what-this-does)
2. [What you need before you start](#2-what-you-need-before-you-start)
3. [Command line survival kit](#3-command-line-survival-kit)
4. [One-time setup](#4-one-time-setup)
5. [Prepare your three input files](#5-prepare-your-three-input-files)
6. [Practice run with the example data](#6-practice-run-with-the-example-data)
7. [Dock your own protein and ligands](#7-dock-your-own-protein-and-ligands)
8. [Follow your job](#8-follow-your-job)
9. [Look at the results](#9-look-at-the-results)
10. [Something went wrong](#10-something-went-wrong)
11. [Glossary](#11-glossary)

---

## 1. What this does

Normal docking lets a ligand move freely in the binding site. **Anchored docking** is for a series of related ligands that share a common core (for example, a set of analogues from one chemical series) when you already have a **template**: a structure of one ligand in the protein, such as a crystal structure.

The workflow finds the part of each new ligand that matches the template, picks a few atoms of that shared core as **anchors**, and gently pulls those atoms toward where the template has them. The rest of the molecule is free to find its best fit.

You provide three files:

| File | What it is |
|---|---|
| **Receptor** (`.pdb`) | The protein. |
| **Template ligand** (`.sdf`) | A ligand already placed correctly in the protein. |
| **Input ligands** (`.sdf`) | The ligand(s) you want to dock. One or many. |

You get back, for every input ligand, a folder containing the docked pose (`ligand_docked.sdf`) and a script to view it in PyMOL.

---

## 2. What you need before you start

- [ ] An account on the HPC and the ability to log in to a **login node** with a terminal (ask your group or IT if you are unsure how).
- [ ] A way to copy files between your computer and the HPC (for example, the file-transfer tool your institute recommends, or `scp`; see [section 5.4](#54-copy-the-files-to-the-hpc)).
- [ ] Your three input files (see [section 5](#5-prepare-your-three-input-files)).
- [ ] [PyMOL](https://pymol.org/) on your own computer to look at results (optional but recommended).

---

## 3. Command line survival kit

The **terminal** is a window where you type instructions and press **Enter** to run them. The text in a grey box in this guide is what you type. Do **not** type the `$` at the start of a line if you see one; this guide leaves it out.

Six things to know:

| You type | What it does |
|---|---|
| `pwd` | Shows which folder you are in. |
| `ls` | Lists the files in the current folder. |
| `cd folder_name` | Moves into a folder. |
| `cd ..` | Moves up one folder. |
| `mkdir folder_name` | Creates a new folder. |
| `cat file_name` | Prints the contents of a text file. |

Tips that save a lot of time:

- **Press Tab** while typing a file or folder name and the terminal completes it for you. Press Tab twice to see the options.
- **Press the Up arrow** to bring back the previous command.
- **Case matters.** `Protein.pdb` and `protein.pdb` are different files.
- **Avoid spaces in file names.** Use `my_protein.pdb`, not `my protein.pdb`.
- **Pasting:** in most terminals, paste with Ctrl+Shift+V (Linux/Windows) or Cmd+V (Mac). Ctrl+C does not copy in a terminal; it **cancels** the running command.
- If a command prints nothing and returns to the prompt, that usually means it worked.

---

## 4. One-time setup

Do this section **once**. Next time, skip to [section 7](#7-dock-your-own-protein-and-ligands).

### 4.1 Log in to the HPC

Open a terminal on your computer and log in to the login node the way you normally do. When you see a prompt that looks like `[yourname@login006 ~]$`, you are on the HPC.

### 4.2 Check that pixi is installed

Pixi is the tool that installs all the docking software for you. Type:

```bash
pixi --version
```

- If you see something like `pixi 0.83.0`, continue to 4.3.
- If you see `command not found`, install it by typing:

  ```bash
  curl -fsSL https://pixi.sh/install.sh | bash
  ```

  Then **close the terminal, log in again**, and re-run `pixi --version`.

The docking assistant finds pixi automatically, so there is nothing else to configure.

### 4.3 Download this workflow

Choose a place for your work (your home or project folder, **not** someone else's shared folder), then download the workflow:

```bash
cd ~
git clone https://github.com/Chemical-Biology-STP/anchored_docking_workflow.git
cd anchored_docking_workflow
```

`cd ~` goes to your home folder. After the last command, your prompt shows you are inside `anchored_docking_workflow`. You can check with `pwd`.

### 4.4 Install the software

Still inside `anchored_docking_workflow`:

```bash
pixi install --all
```

This downloads and installs RDKit, Meeko, AutoGrid and ADFRsuite. It takes a few minutes and prints a lot of text. Wait until you get your prompt back and see a line such as `The following environments have been installed: default, adfr`.

You do not need to install AutoDock-GPU: the HPC provides it as a **module** that the job script loads for you.

Setup is finished.

---

## 5. Prepare your three input files

This is the step where most problems start, so go through it carefully. Put all three files in a folder **on your own computer** first, with simple names, for example:

```
receptor.pdb
template.sdf
ligands.sdf
```

### 5.1 The receptor (`receptor.pdb`)

Start from the same protein structure that your template ligand came from (for example, the PDB entry of the crystal structure). Then check:

- **Only the protein.** Remove the ligand, waters, ions and crystallisation additives unless you deliberately want to keep them. In PyMOL, you can do this with:

  ```
  remove solvent
  remove organic
  ```

- **Hydrogen atoms must be present.** The receptor is processed *without* adding hydrogens, so a receptor with no hydrogens will give poor results. Many crystal structures from the PDB have none. To add them in PyMOL:

  ```
  h_add
  ```

  PyMOL's hydrogen placement is quick but basic (for example, it makes no careful choice of histidine protonation states). If your project needs more care, protonate the protein with your preferred tool first.
- **One structure only.** If the file holds several models or alternate conformations, keep one.
- Save it in PDB format: in PyMOL, `File > Export Molecule...`, choose `PDB`.

### 5.2 The template ligand (`template.sdf`)

The template is a ligand that is **already in the right place in the receptor**. Typically this is the ligand from the crystal structure.

- Its coordinates must be in **the same coordinate system as your receptor.** The easiest way to guarantee this is to export both receptor and ligand from the **same** PDB file, without moving or aligning either of them afterwards.
- In PyMOL: select the ligand (for example `select lig, organic`), then `File > Export Molecule...`, choose `lig` and the format `SDF`.
- Only the **first** molecule in the file is used as the template.
- The template must contain the common core that your input ligands share.

### 5.3 The input ligands (`ligands.sdf`)

- An SDF file with one or more molecules. Each molecule is docked separately, in file order.
- Hydrogens are not required; they are added automatically.
- The molecules must be valid chemistry (RDKit has to be able to read them). If you can open the file in a chemistry program, it is usually fine.
- **Each input ligand needs to share a substructure with the template.** The workflow looks for that shared core automatically. If a ligand shares very little with the template, it cannot be docked this way (see [section 10](#10-something-went-wrong)).
- Give each molecule a distinct name in the first line of its SDF block. This makes the results easy to tell apart.

### 5.4 Copy the files to the HPC

Put the files in a new folder inside the workflow folder. On the HPC terminal:

```bash
cd ~/anchored_docking_workflow
mkdir my_project
```

Now copy your three files from your computer into `~/anchored_docking_workflow/my_project/`. Two ways:

- **With a file-transfer program** (graphical, drag and drop): connect to the HPC and drag the files into the `my_project` folder.
- **With `scp`**: open a *new* terminal on **your own computer** (not the HPC) and type, replacing the parts in capitals:

  ```bash
  scp receptor.pdb template.sdf ligands.sdf YOURNAME@LOGIN_NODE_ADDRESS:~/anchored_docking_workflow/my_project/
  ```

Back on the HPC, check that they arrived:

```bash
ls my_project
```

You should see `ligands.sdf  receptor.pdb  template.sdf`.

---

## 6. Practice run with the example data

Before using your own data, run the small example that ships with the workflow. If this works, your setup is correct.

All commands are typed on the HPC, inside the workflow folder:

```bash
cd ~/anchored_docking_workflow
pixi run dock
```

This starts the **interactive docking assistant**. It asks one question at a time. Each question shows a suggested answer in `[square brackets]`: **press Enter to accept it**. Press **Tab** to complete file names, and **Ctrl+C** at any time to quit (nothing is run until you confirm at the end).

Answer like this:

| Step | Question | What to type |
|---|---|---|
| 1 | Receptor PDB file | `test/Bace1_protein.pdb` |
| 2 | Template ligand SDF file | `test/crystal_ligand.sdf` |
| 3 | Ligand SDF file | `test/CAT-13f.sdf` |
| 4 | Use the recommended docking settings? | Enter |
| 4 | Delete intermediate files? | Enter |
| 5 | Choose a number (how to run) | Enter (`1`, submit to the cluster) |
| 5 | Output folder | `practice_output` |
| 5 | Job name, Time limit, Change partition...? | Enter for each |
| 6 | Submit the job to the queue? | Enter |

After each of the first three answers the assistant checks your file and prints lines starting with `[ok]`, `[!]` (a warning) or `[x]` (a problem). For the example you should only see `[ok]` lines.

On this cluster you may then see a submission estimate and a question `Submit now? [Y/n]`. Type `Y` and press Enter. The assistant reports `Submitted job 12345678`, then asks `Watch the job until it finishes?`. Press Enter to say yes.

It prints a line whenever the job changes state (waiting in the queue, running) and finishes with a table of results. The example takes only a few minutes once it starts. When it is done, the pose is in `practice_output/0/ligand_docked.sdf`.

If you press **Ctrl+C** while watching, the assistant stops watching but **your job keeps running**. See [section 8](#8-follow-your-job) to check on it later.

---

## 7. Dock your own protein and ligands

Once the practice run works:

```bash
cd ~/anchored_docking_workflow
pixi run dock
```

Always start from this folder (`anchored_docking_workflow`). Paths you type can be absolute (`/home/...`) or relative to this folder, for example `my_project/receptor.pdb`.

### Steps 1 to 3: your three files

Type the path to each file (Tab completes it). After each one the assistant checks it and explains any problem:

| Check | What it tells you |
|---|---|
| Receptor has hydrogens | `[!]` if there are none. Fix this before continuing (see [5.1](#51-the-receptor-receptorpdb)); results will be poor otherwise. |
| Receptor contains other molecules / several models | `[!]` as a reminder to clean the file. |
| Template has 3D coordinates | `[x]` if not; you must give a ligand that is positioned in the receptor. |
| Template sits in the receptor | `[!]` if most of its atoms are far from the protein. This usually means the template and receptor are in different coordinate systems. Press Ctrl+C and re-export both from the same structure ([5.2](#52-the-template-ligand-templatesdf)). |
| Template file holds one molecule | `[!]` if there are several; only the first is used. |
| Each ligand shares a core with the template | The assistant runs the real core and anchor search on every ligand. Ligands with no shared core, or only a tiny one, are listed. |
| Ligand is a single molecule | `[!]` if a ligand has several fragments, for example a salt with a counter-ion. Remove the extra fragments in your chemistry program. |
| Ligand atoms are supported | `[!]` if a ligand contains an element AutoDock has no atom type for (such as selenium). |

If some ligands are flagged, the assistant offers to **leave them out and dock the rest**. This is usually what you want: a single ligand that cannot be docked would otherwise stop the whole job when it is reached. The remaining ligands are saved in a new file called `<job name>_ligands.sdf` next to the job script, and the **result folders are numbered within that reduced list**. The review screen shows which folder number is which ligand.

### Step 4: docking settings

You are first offered the recommended settings: 3 anchors (or fewer if the weakest ligand has fewer usable ones), random anchor choice, and a 20 Å search box. Press Enter to accept them. If you answer `no`, you are asked for:

| Question | What it means | Advice |
|---|---|---|
| Anchor selection | `random` picks anchors at random from the shared core; `area` tries many random sets and keeps the three that span the largest triangle | Keep `random` to start. Try `area` if poses are rotated about the core. `area` always uses exactly 3 anchors. |
| Number of anchor atoms (`random` only) | How many shared-core atoms are pulled towards the template | Keep `3`. More anchors hold the ligand more rigidly. The most the software supports is 17, and the assistant refuses larger numbers. |
| Number of random attempts (`area` only) | How many anchor sets to try | Keep `100`. |
| Search box size | Size (Å) of the cube around the template in which docking happens | Keep `20`. The assistant warns if the box is not larger than the template, and refuses more than 95 Å. Run time and disk use grow with the cube of the size, so a 40 Å box needs about 8 times the disk space of a 20 Å box. |
| Random seed | A number that makes the anchor choice repeatable | Enter a whole number of 0 or higher, e.g. `1`, to be able to rerun and get the same anchors. Leave blank otherwise. |
| Expert options | A manual core pattern (SMARTS) and debug output | Say `no` unless you know what a SMARTS pattern is. |

Finally you are asked whether to **delete intermediate files** after each ligand. They take about 30 MB per ligand. The suggestion is `yes` for 10 or more ligands and `no` for fewer. Choose `no` if you want to look at scores later, because the score file (`ligand.dlg`) is deleted with the rest.

### Step 5: how to run it

| Choice | When to use it |
|---|---|
| **1. Submit to the cluster queue** | Normal choice. The assistant writes a job script, submits it, and can watch it. |
| **2. Run now in this terminal** | Only if you are already on a GPU machine with the AutoDock-GPU module loaded. It checks for both and tells you if something is missing. |
| **3. Only save the job script** | If you want to look at, or edit, the script before submitting it yourself with `sbatch`. |

Then you choose the **output folder**, where results are written (one numbered sub-folder per ligand). The assistant protects you here: it refuses a folder that contains your input files, your home folder or this workflow, and it asks before deleting an existing non-empty folder. Avoid spaces in the name.

For cluster jobs you are also asked for:

- **Job name:** anything recognisable.
- **Time limit:** the suggestion allows 10 minutes per ligand (at least 30 minutes). The example ligand needed about 2 minutes, but bigger proteins and ligands take longer. If the time runs out, the job is stopped, so raise it for large jobs.
- **Change partition, GPUs, modules or pixi location?** Say `no`. The defaults are the GPU partition `ga100`, 1 GPU, modules `GCC/13.2.0` and `AutoDock-GPU/1.5.3-CUDA` (both are required: AutoDock-GPU only works with GCC 13.2), and the pixi folder found automatically.

### Step 6: review

The assistant prints every setting and any warnings, then asks for a final confirmation. **Nothing is submitted until you say yes here.** Read the warnings: they list every `[!]` raised so far.

The job script is saved next to your output folder as `submit_<job name>.sh`, and the log as `<job name>_<jobnumber>.log`. You can print the script with `cat` to see exactly what was run, and keep it as a record of the settings you used.

---

## 8. Follow your job

If you let the assistant watch the job, it tells you when the job starts, runs and finishes, then shows a table of which ligands produced a pose. If the job failed, it also prints the last lines of the log.

If you stopped watching (Ctrl+C), closed the terminal, or chose to only save the script, use these commands. Your job first **waits in the queue**, then **runs**.

```bash
squeue -u $USER
```

The `ST` column shows `PD` (pending, waiting for a free GPU) or `R` (running). If your job is not listed, it has finished (successfully or not).

Look at the end of the log with:

```bash
tail my_project/dock_myjob_12345678.log
```

(Replace the file name with your own; Tab completion helps.) To watch it live, use `tail -f` instead of `tail`, and press **Ctrl+C** to stop watching. That does not cancel the job.

To cancel a job:

```bash
scancel 12345678
```

If you only saved the script, submit it with `sbatch my_project/submit_myjob.sh`.

### How do I know it worked?

- The results table lists every ligand as `docked`, with small numbers in the last two columns (see [9.0](#90-reading-the-results-table)).
- The job state is `COMPLETED`. If one or more ligands failed it shows `FAILED` instead, but the other ligands are still docked, so read the table before concluding that nothing worked.
- The log ends with a line like `Docking complete.`
- There is a folder for each ligand in your output directory, named `0`, `1`, `2`, ... in the same order as in your input SDF (or the reduced list, if you left ligands out).
- Inside each folder is a file `ligand_docked.sdf`.

Check:

```bash
ls my_project/output
ls my_project/output/0
```

---

## 9. Look at the results

Each ligand has its own folder, `output/0`, `output/1`, and so on. The important files are:

| File | What it is |
|---|---|
| `ligand_docked.sdf` | The docked pose. Open it in any molecule viewer. |
| `template.sdf` | A copy of your template ligand. |
| `receptor.pdb` | A copy of your receptor. |
| `view_results.pml` | A script that loads all three in PyMOL with sensible colours. |
| `ligand.dlg` | The full AutoDock-GPU log, with scores. |
| `anchor_positions.json` | The coordinates of the anchor points that were used. |

If you chose to delete intermediate files, only the first four (plus `docking_results.pse` once you have run the PyMOL script) are kept.

### 9.0 Reading the results table

When the job ends, the assistant prints a table with one row per ligand:

| Column | Meaning |
|---|---|
| Folder | The result folder number (`0`, `1`, ...). |
| Result | `docked`, `CHECK POSE` or `FAILED`. |
| Anchor dev (A) | How far the most displaced anchor atom ended up from the position the template asked for. |
| Core RMSD (A) | How far the shared core is from the template's core, averaged over all its atoms. |

The same information is saved in `summary.tsv` in your output folder.

- **`docked` with small numbers** (roughly under 1.5 Å for anchors and under 1 Å for the core in the example) means the ligand sits where the anchoring asked.
- **`CHECK POSE`** means the pose was written but did not satisfy the anchoring: an anchor is more than 2.0 Å from its target, or the core is more than 2.5 Å RMSD from the template. Docking only *pulls* the anchors towards the template, so this can happen. Open the pose in PyMOL before trusting it. These limits are judgements based on a few runs, so look at the numbers as well as the label.
- **`FAILED`** means that ligand could not be docked. The reason is printed in the table and saved in `error.txt` inside its folder. The other ligands are still docked.

To redo only the failed ligands after fixing the cause, run `pixi run dock` again with the **same files and settings** and the **same output folder**. The assistant notices the earlier run and asks `Resume it?`. Answer yes: finished ligands are skipped. If you changed the files or settings, resuming is refused, and you can choose a new output folder instead.

### 9.1 Copy the results to your computer

Copy the whole output folder. In a terminal on **your own computer**:

```bash
scp -r YOURNAME@LOGIN_NODE_ADDRESS:~/anchored_docking_workflow/my_project/output .
```

(`-r` means "the whole folder"; the final `.` means "into the folder I am currently in".) Or use your file-transfer program.

### 9.2 Open in PyMOL

1. Start PyMOL.
2. In the PyMOL command line (at the bottom of the window), move into one ligand's folder and run the script. Use the real path on your computer:

   ```
   cd /path/to/output/0
   @view_results.pml
   ```

3. You see the receptor (grey cartoon), the template ligand (cyan) and the docked ligand (magenta).

What to look for:

- The docked ligand's **core should sit on top of the template's core**. That is the point of anchoring.
- The rest of the ligand should sit in the pocket without clashing with the protein.
- If the pose looks odd (flipped, far from the template, or clashing), try again with a different anchor strategy or seed (see [section 10](#10-something-went-wrong)).

### 9.3 Scores

Open `ligand.dlg` and look for the ranking lines:

```bash
grep RANKING my_project/output/0/ligand.dlg
```

Each line has several columns; the fourth is the estimated binding energy in kcal/mol (more negative is better). Keep in mind that the anchors add an artificial bias to the grid maps, so these numbers are only a rough guide and should not be compared with scores from unbiased docking. Use your eyes on the pose first.

---

## 10. Something went wrong

Most problems with input files are caught by the assistant before anything is submitted, so read its `[!]` and `[x]` lines carefully. If a job fails, look at the **end** of the log file (`tail <jobname>_<jobnumber>.log`) and find the message in this table.

| Message or symptom | What it means and what to do |
|---|---|
| `This tool must be started through pixi` | You started the script directly. Run `pixi run dock` from the `anchored_docking_workflow` folder. |
| `ADFRsuite ... is not installed yet` | You ran `pixi install` without `--all`. Run `pixi install --all` once, then `pixi run dock` again. |
| `You do not have permission to write in ...` | Choose an output folder inside your own home or project folder. |
| `This run may need about N GB but only M GB is free` | Free up space, choose a folder on a larger disk, delete intermediate files, or use a smaller box. |
| `Could not reach SLURM to check the job` | The cluster scheduler did not answer. Your job is not affected. Check on it with `squeue -u $USER` in a minute. |
| `pixi: command not found` | pixi is not on your path. Run `pixi --version`; if that fails, install pixi ([4.2](#42-check-that-pixi-is-installed)). Then run `pixi run dock` again so the assistant can find it. |
| `AutoDock-GPU binary not found ... PATH` | The modules were not loaded. For cluster jobs, check with `cat submit_<jobname>.sh` that it contains the two `module load` lines; if not, run `pixi run dock` again and keep the default modules. For "Run now", run `module load GCC/13.2.0 AutoDock-GPU/1.5.3-CUDA` first. |
| `GLIBCXX_3.4.29 not found` | `GCC/13.2.0` was not loaded before AutoDock-GPU. Make sure both modules are listed, with `GCC/13.2.0` first. |
| `No GPU is visible on this machine` | You chose "Run now" on a machine without a GPU, such as the login node. Choose "Submit to the cluster queue" instead. |
| `[x]` "could not generate 3D coordinates" | RDKit could not build a 3D shape for this ligand, usually because of a strained ring system or impossible stereochemistry. Check the structure in your chemistry program, or let the assistant leave it out. |
| Result `CHECK POSE` | The pose did not satisfy the anchoring. See [9.0](#90-reading-the-results-table). |
| Result `FAILED` for a ligand | Read the reason in the table and the `error.txt` in that ligand's folder. Fix the cause, then run `pixi run dock` again with the same output folder and answer yes to resume. |
| `Cannot resume: these inputs or settings differ` | Resuming only works when the files and docking settings are identical to the earlier run. Use a new output folder, or delete the old one. |
| `[x]` "no shared core with the template" | The ligand has nothing in common with the template, so it cannot be anchored. Check it belongs to the same series, or let the assistant leave it out. |
| `[!]` ligand "shares only N heavy atoms" | The ligand is very different from the template. Anchors on such a small core are probably meaningless. |
| `Not enough anchor points found.` | The shared core of a ligand has fewer usable anchor atoms than you asked for. Lower the number of anchors. The assistant normally prevents this by suggesting a safe number. |
| `Failed to load molecule files`, or an error mentioning `NoneType` | RDKit could not read a molecule in one of the SDF files. Open the file in a chemistry program, find the problem molecule and fix or remove it. |
| `File not found` | A typo in the path, or you are not in the `anchored_docking_workflow` folder. Check with `pwd` and `ls`. |
| Job stays `PD` for a long time | The GPU partition is busy. Wait, or ask for less time (shorter jobs are often scheduled sooner). |
| Job disappeared and there are only some output folders | The time limit ran out. Increase the time limit and rerun (or put the remaining ligands in a new SDF file). |
| Docked ligand is far from the template | The template and receptor are not in the same coordinate system. Re-export both from the **same** structure (see [5.2](#52-the-template-ligand-templatesdf)). |
| Docked ligand clashes with the protein | The receptor may lack hydrogens or contain leftover molecules. Re-check [5.1](#51-the-receptor-receptorpdb). A larger search box or a different seed may also help. |
| Results differ between runs | `random` anchor mode picks different anchors each run. Set a random seed to repeat a run exactly. |

If you are still stuck, send your colleague or the maintainers: the **last 30 lines of the log** (`tail -n 30 <jobname>_<jobnumber>.log`) and the **contents of the job script** (`cat submit_<jobname>.sh`). Those two items answer most questions.

---

## 11. Glossary

- **Terminal / command line:** a text window for giving instructions to the computer.
- **HPC:** high-performance computing cluster, a shared set of powerful computers.
- **Login node:** the computer you connect to first. It is for preparing and submitting work, not for heavy calculations.
- **SLURM:** the system that queues jobs and hands them out to compute nodes. `sbatch` submits a job, `squeue` shows the queue, `scancel` cancels a job.
- **GPU:** a graphics processor. AutoDock-GPU uses it to dock quickly.
- **Module:** a pre-installed program on the HPC that you switch on with `module load`.
- **pixi:** a package manager that installs the Python and docking tools into a project folder, so you do not install anything by hand.
- **Receptor:** the protein you dock into.
- **Template ligand:** a ligand already correctly placed in the receptor.
- **Anchor:** an atom of the shared core that is pulled toward the template's position for that atom.
- **SDF / PDB:** common file formats for molecules (ligands / proteins).
- **SMARTS:** a text notation for describing a molecular substructure (advanced).
