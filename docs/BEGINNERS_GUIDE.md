# Beginner's Guide: Anchored Docking on the HPC

This guide takes you from "I have a protein and some ligands" to "I have docked poses I can look at", assuming you have never used a command line before. Follow the steps in order and copy the commands exactly.

> **Status:** this workflow is under active development and has not been experimentally validated. Treat results as hypotheses to check, not as facts.

## Contents

1. [What this does](#1-what-this-does)
2. [What you need before you start](#2-what-you-need-before-you-start)
3. [Command line survival kit](#3-command-line-survival-kit)
4. [One-time setup](#4-one-time-setup)
5. [Prepare your input files](#5-prepare-your-input-files)
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

You provide:

| File | What it is |
|---|---|
| **Structure of the complex** (`.cif` or `.pdb`) | The protein with its bound ligand, for example downloaded from the PDB. The assistant extracts the template ligand and the cleaned protein from it. |
| **Input ligands** (`.sdf`) | The ligand(s) you want to dock. One or many. |

If you already have a cleaned protein and a template ligand as separate files, you can give those instead of the complex: a **receptor** (`.pdb`) and a **template ligand** (`.sdf`) already placed correctly in the protein.

You get back, for every input ligand, a folder containing the docked pose (`ligand_docked.sdf`) and a script to view it in PyMOL.

---

## 2. What you need before you start

- [ ] An account on the HPC and the ability to log in to a **login node** with a terminal (ask your group or IT if you are unsure how).
- [ ] A way to copy files between your computer and the HPC (for example, the file-transfer tool your institute recommends, or `scp`; see [section 5.4](#54-copy-the-files-to-the-hpc)).
- [ ] Your input files: a structure of the protein with its ligand, and the ligands to dock (see [section 5](#5-prepare-your-input-files)).
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

## 5. Prepare your input files

There are two ways to provide the protein and the template ligand. **Most people should use the first**: you give the assistant the structure file of the protein-ligand complex, and it does the clean-up for you. Your ligands to dock (section 5.3) are needed either way.

Put your files in a folder **on your own computer** first, with simple names, for example:

```
complex.cif        (route A: the structure of the protein with its ligand)
ligands.sdf        (the ligands you want to dock)
```

or, for route B:

```
receptor.pdb
template.sdf
ligands.sdf
```

### 5.0 Easiest: start from the complex (route A)

Download the structure of your protein with its bound ligand from the [PDB](https://www.rcsb.org): open the entry (for example `1M17`), choose **Download Files**, then **mmCIF Format** (a `.cif` file). A `.pdb` file also works.

Nothing needs to be cleaned by hand. The assistant reads the file (your copy is never changed) and:

- finds the molecules in it and **asks which one is the template ligand**, suggesting the most likely one;
- **removes waters and crystallisation additives** (sulfate, glycerol, PEG and similar) automatically, and tells you what it removed;
- **asks which protein chains to keep** if there are several, suggesting the ones next to the template;
- **asks about each ion, cofactor or other molecule** near the template, so you decide whether it belongs to the receptor (for example a metal in the active site);
- keeps only the first model and the first alternate conformation of each atom;
- writes the cleaned protein and the template ligand into your output area for the job;
- **prepares the protein**: it reports chain gaps and missing atoms, and if the structure has no hydrogens (most structures have none) it protonates the protein at a pH you choose (default 7.4); see [Protein preparation](#protein-preparation).

Things to know:

- The template is taken from the coordinates only; the structure file does not say which bonds are double, so the template is drawn with single bonds in PyMOL. This does not change the anchors or the docking (checked against a template with correct bonds).
- Ions and cofactors that you keep may lack the partial charges AutoDock uses, which weakens their effect.
- Protonation is predicted (PDB2PQR with PROPKA), not measured. The assistant lists the residues next to the template that deserve a look, and the full record is saved in `<job name>_protonation_report.txt`.

### Protein preparation

Docking needs hydrogens on the protein, and which hydrogens matter. A histidine can carry a hydrogen on either ring nitrogen or on both; an aspartate in a buried pocket may be neutral instead of negative. The assistant does this for you in step 2, for both routes:

1. **Gaps and missing atoms are reported first.** A *gap* means residues are missing from the chain in the structure file (common in loops that the crystal did not resolve). A gap far from the template is only noted; one within 10 Å of the template, or a residue next to the template that lacks atoms, is shown as a warning, because the pocket shape may be wrong there. The assistant cannot invent missing loops; if the gap matters, model it elsewhere first.
2. **If the protein has no hydrogens**, it asks whether to protonate it and at which pH. The default is **7.4**. It uses PDB2PQR with PROPKA: PROPKA predicts the pKa of every titratable residue from its surroundings, and PDB2PQR then chooses each histidine's form, protonates aspartates, glutamates and the chain ends accordingly, and places the hydrogens so that hydrogen bonds fit. It also rebuilds missing side-chain atoms and may turn the amide of an asparagine or glutamine, or a histidine ring, by 180° when that gives better hydrogen bonds. This takes under a minute for a typical protein.
3. **What you see afterwards:** the residues within 6 Å of the template whose state is unusual or whose pKa is close to the pH (these are the ones worth a second look), and any residue PDB2PQR could not handle cleanly. The full record is saved as `<job name>_protonation_report.txt`: all unusual states, flipped side chains, rebuilt atoms, gaps, and a pKa for every titratable residue near the template. The protonated protein is saved as `<job name>_receptor.pdb` and is what is docked.

If the protein **already has hydrogens**, they are used as they are and not checked. If you answer **no** to protonation, the assistant offers ADFR's quicker rule-based hydrogens, which treat every histidine as positively charged and give gap ends artificial charges. Ions and cofactors that you keep are not given hydrogens.

If PDB2PQR fails for some reason, the assistant says so and falls back to the ADFR offer.

### 5.1 The receptor (`receptor.pdb`)

Route B only. Use this if you already have a cleaned protein, or you want full control. Start from the same protein structure that your template ligand came from (for example, the PDB entry of the crystal structure). Then check:

- **Only the protein.** Remove the ligand, waters, ions and crystallisation additives unless you deliberately want to keep them. In PyMOL, you can do this with:

  ```
  remove solvent
  remove organic
  ```

- **Hydrogen atoms must be present**, or let the assistant add them. A receptor with no hydrogens gives poor results, and many crystal structures from the PDB have none. The assistant offers to protonate the protein for you (PDB2PQR with PROPKA, at pH 7.4 unless you choose otherwise). **A receptor that already has hydrogens is used exactly as it is**, including its histidine forms and charged or neutral acids; they are not checked. To add hydrogens yourself in PyMOL:

  ```
  h_add
  ```

  PyMOL's hydrogen placement is quick but basic (for example, it makes no careful choice of histidine protonation states), so letting the assistant protonate a file without hydrogens is the better route.
- **One structure only.** If the file holds several models or alternate conformations, keep one.
- Save it in PDB format: in PyMOL, `File > Export Molecule...`, choose `PDB`.

### 5.2 The template ligand (`template.sdf`)

Route B only.

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

Now copy your files from your computer into `~/anchored_docking_workflow/my_project/`. Two ways:

- **With a file-transfer program** (graphical, drag and drop): connect to the HPC and drag the files into the `my_project` folder.
- **With `scp`**: open a *new* terminal on **your own computer** (not the HPC) and type, replacing the parts in capitals:

  ```bash
  scp complex.cif ligands.sdf YOURNAME@LOGIN_NODE_ADDRESS:~/anchored_docking_workflow/my_project/
  ```

Back on the HPC, check that they arrived:

```bash
ls my_project
```

You should see your files, for example `ligands.sdf  complex.cif` (route A) or `ligands.sdf  receptor.pdb  template.sdf` (route B).

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
| 1 | What do you have? (Choose a number) | `2` (the example comes as a separate protein and template, see [5.1](#51-the-receptor-receptorpdb)) |
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

### Steps 1 to 3: your protein, template and ligands

The first question is **What do you have?** Answer `1` if you have a structure file of the protein with its ligand (route A, recommended), or `2` for a cleaned protein plus a separate template file (route B).

**Route A (answer 1).** Type the path to your `.cif` or `.pdb` file. The assistant then asks, in this order:

| Question | What to do |
|---|---|
| Template number | A table lists the molecules in the file with how many protein atoms each touches. The most likely template is suggested. Press Enter to accept it, type another number, or `q` to choose another file. If nothing in the file is ligand-sized, it tells you and asks for another file. |
| Chains to keep | Only asked if there are several protein chains. The chains next to the template are suggested; press Enter to accept, or type letters separated by spaces. |
| Keep *name* (ion/cofactor, N atoms, X Å from the template)? | One question for each ion, cofactor or other molecule within 12 Å of the template. Molecules within 5 Å of the template default to *keep*; others to *remove*. Say yes for something that belongs to the receptor, such as a catalytic metal. Molecules farther than 12 Å, waters, additives and sugars are removed without asking, and listed afterwards. If there are many, similar molecules are asked about together. |
| Protonate the protein with PDB2PQR and PROPKA? / pH | Only asked if the structure has no hydrogens. Answer yes (the default) and accept pH 7.4 unless your target works at another pH (for example an acidic compartment). If you say no, you are offered ADFR's quicker rule-based hydrogens instead. |

The assistant then prints **what was done with the structure**: which template and chains were used and everything that was removed. Read this list. The cleaned protein and the template are saved next to your results as `<job name>_receptor.pdb` and `<job name>_template.sdf`, and you can open them in PyMOL to check.

**Route B (answer 2).** Type the path to the protein PDB file, then to the template SDF file. If the protein has no hydrogens you are asked about protonation (see [Protein preparation](#protein-preparation)); your own file is never changed, and the protonated copy is saved with the results.

**Both routes, then:** type the path to your ligands. After each file the assistant checks it and explains any problem:

| Check | What it tells you |
|---|---|
| Receptor has hydrogens | `[!]` if there are none, followed by the offer to protonate. If you decline both offers, results will be poor. |
| Receptor contains other molecules / several models (route B) | `[!]` as a reminder to clean the file. Route A removes these itself. |
| Template has 3D coordinates | `[x]` if not; you must give a ligand that is positioned in the receptor. |
| Template sits in the receptor | `[!]` if most of its atoms are far from the protein. This usually means the template and receptor are in different coordinate systems. Press Ctrl+C and re-export both from the same structure ([5.2](#52-the-template-ligand-templatesdf)). |
| Template file holds one molecule | `[!]` if there are several; only the first is used. |
| Each ligand shares a core with the template | The assistant runs the real core and anchor search on every ligand. Ligands with no shared core, or only a tiny one, are listed. |
| Ligand is a single molecule | `[!]` if a ligand has several fragments, for example a salt with a counter-ion. Remove the extra fragments in your chemistry program. |
| Ligand atoms are supported | `[!]` if a ligand contains an element AutoDock has no atom type for (such as selenium). |

If some ligands are flagged, the assistant offers to **leave them out and dock the rest**. This is usually what you want: a single ligand that cannot be docked would otherwise stop the whole job when it is reached. The remaining ligands are saved in a new file called `<job name>_ligands.sdf` next to the job script, and the **result folders are numbered within that reduced list**. The review screen shows which folder number is which ligand.

### Step 4: docking settings

**What an anchor is.** For each ligand the assistant finds the atoms it shares with the template (the maximum common substructure) and uses some of them as *anchors*. Each anchor atom is gently pulled towards the position of its matching atom in the template. The parts of the ligand that differ from the template are free.

You are first offered the recommended settings: **every usable atom of the shared core as an anchor, chosen to be spread across the core** (so the whole core is held in place), and a 20 Å search box. Press Enter to accept them. Two limits apply. An atom with a symmetric twin in the core is never used as an anchor. And AutoDock-GPU supports 32 atom types per ligand in all, counting the anchors and the ordinary types of the ligand's other atoms, so a ligand with a very large core gets a few anchors fewer; the results table says when that happens (for example `only 29 of 30 anchors could be used`). If you answer `no`, you are asked for:

| Question | What it means | Advice |
|---|---|---|
| Anchor selection | `spread` takes up to a maximum number of shared-core atoms that are furthest apart from each other (no randomness); `random` picks a number of anchors at random; `area` tries many random sets and keeps the three that span the largest triangle | Keep `spread`. Use `random` or `area` only to compare with earlier runs. `area` always uses exactly 3 anchors. |
| Maximum number of anchor atoms (`spread` only) | The most anchors any ligand gets; a ligand with a smaller usable core uses all it has | Keep the maximum. Fewer anchors let the core move more. |
| Number of anchor atoms (`random` only) | How many shared-core atoms are pulled towards the template | Try `3` to start. More anchors hold the ligand more rigidly. The assistant refuses numbers that cannot fit AutoDock-GPU's 32 atom types; a ligand that needs more fails with an explanation. |
| Number of random attempts (`area` only) | How many anchor sets to try | Keep `100`. |
| Search box size | Size (Å) of the cube around the template in which docking happens | Keep `20`. The assistant warns if the box is not larger than the template, and refuses more than 95 Å. Run time and disk use grow with the cube of the size, so a 40 Å box needs about 8 times the disk space of a 20 Å box. |
| Random seed (not asked in `spread` mode, which has no random element) | A number that makes the anchor choice repeatable | Enter a whole number of 0 or higher, e.g. `1`, to be able to rerun and get the same anchors. Leave blank otherwise. |
| Expert options | A manual core pattern (SMARTS) and debug output | Say `no` unless you know what a SMARTS pattern is. |

Finally you are asked whether to **delete intermediate files** after each ligand. They take about 30 MB per ligand. The suggestion is `yes` for 10 or more ligands and `no` for fewer. Deleting them is safe for scores: they are kept in `results.csv` and in the SDF files.

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
- **Change partition, GPUs, modules or pixi location?** Say `no`. The defaults are the GPU partition `ga100`, 1 GPU (for the H100 nodes enter partition `gh100`; the assistant then suggests the reservation `h100`, which that partition requires and which only some accounts may use. The assistant asks the cluster whether your account can use the partition you chose, before it writes the job script, and tells you if you cannot), modules `GCC/13.2.0` and `AutoDock-GPU/1.5.3-CUDA` (both are required: AutoDock-GPU only works with GCC 13.2), and the pixi folder found automatically.

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

The job does not write a PyMOL session. It writes PyMOL **scripts** (text files of PyMOL commands); the session file appears when you run a script in PyMOL.

In your output folder you will find:

| File | What it is |
|---|---|
| `results.csv` | **One row per ligand**: status, binding energy (and its parts), anchor deviation, core RMSD, failure reason. Opens in Excel. |
| `top_poses.sdf` | The best pose of every docked ligand in one file, with the scores stored as data fields. |
| `all_poses.sdf` | Every pose the docking kept (one representative per cluster of similar poses), for every ligand. `pose_rank` 1 is the best. |
| `view_all.pml` | A PyMOL script that opens the protein, the template and **all** docked ligands in one session. |
| `receptor_with_hydrogens.pdb`, `template_with_hydrogens.sdf` | The protein and template that the PyMOL scripts load, with every hydrogen. |
| `summary.tsv`, `run_settings.json` | The per-ligand status, and the settings used (needed to resume). |
| `0/`, `1/`, ... | One folder per ligand (see below). |

Inside each ligand's folder:

| File | What it is |
|---|---|
| `ligand_docked.sdf` | The docked pose or poses for this ligand, as complete molecules with all hydrogens. Open it in any molecule viewer. |
| `view_results.pml` | A PyMOL script that opens just this ligand with the protein and template. |
| `pose_check.json`, `anchor_positions.json` | The pose check numbers, and the coordinates of the anchor points that were used. |
| `error.txt` | Only if this ligand failed. |
| `ligand.dlg` and other files | The raw AutoDock-GPU output and intermediate files. Deleted if you chose to delete intermediate files. |

### Hydrogens

AutoDock works with the non-polar hydrogens merged into their carbon atoms, so its own files are not what you want to look at. Everything you open in PyMOL has **all hydrogens back**: the docked ligands, the protein and the template. If your protein file had no hydrogens and you let the assistant add them, the protein you see is the one that was docked, plus its non-polar hydrogens. The polar hydrogens, which matter for hydrogen bonds, are exactly where they were during docking.

### 9.0 Reading the results table

When the job ends, the assistant prints a table with one row per ligand:

| Column | Meaning |
|---|---|
| Folder | The result folder number (`0`, `1`, ...). |
| Result | `docked`, `CHECK POSE` or `FAILED`. |
| Energy (kcal/mol) | The docking score of the best pose; more negative is better (see [9.3](#93-scores)). |
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
2. In the PyMOL command line (at the bottom of the window), move into the output folder you copied and run the script. Use the real path on your computer:

   ```
   cd /path/to/output
   @view_all.pml
   ```

   (To look at one ligand only, move into its folder, for example `cd /path/to/output/0`, and run `@view_results.pml`.)

3. You see the protein (grey cartoon, with the residues around the ligands as sticks), the template ligand (cyan) and every docked ligand (magenta; poses the pose check flagged are orange). The scripts save a `docking_all.pse` (or `docking_results.pse`) session file that you can reopen later without the script.

PyMOL must be started from, or moved into, the folder holding the script, because the script refers to the other files by short relative names. If PyMOL says it cannot find a file, check the folder with `pwd` in the PyMOL command line.

What to look for:

- The docked ligand's **core should sit on top of the template's core**. That is the point of anchoring.
- The rest of the ligand should sit in the pocket without clashing with the protein.
- A ligand with several poses holds them as separate *states*: use the arrows at the bottom right of the PyMOL window to step through them.
- If a pose looks odd (flipped, far from the template, or clashing), try again with a different anchor strategy or seed (see [section 10](#10-something-went-wrong)).

### 9.3 Scores

The binding energies are in `results.csv` (column `free_energy`, in kcal/mol; more negative is better), in the results table printed at the end of the job, and in `top_poses.sdf`, where each molecule carries `free_energy`, `intermolecular_energy`, `internal_energy` and `cluster_size`. `cluster_size` is how many of the 20 docking runs ended in that pose.

Keep in mind that the anchors add an artificial bias to the grid maps, so the energies are only a rough guide, are best used to rank ligands **within this run**, and should not be compared with scores from unbiased docking. Use your eyes on the pose first.

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
| `No molecule in this file is big enough to be a ligand` | The structure has no bound ligand (an apo structure), or only ions and tiny molecules. Choose a structure that contains your template ligand. |
| `Could not read this file as a PDB or mmCIF structure` | The file is damaged or is not a structure file. Download it again from rcsb.org (mmCIF format). |
| `No protein chain was found in this structure` | The file holds no protein, for example only DNA or a ligand. Check you downloaded the right entry. |
| The suggested template is the wrong molecule | Type the number of the one you want in the table instead of pressing Enter. The table shows how many protein atoms each molecule touches, which helps to tell a bound drug from a stray molecule. |
| Cofactor or ion you wanted is missing from the receptor | You answered no to it, or it was more than 12 Å from the template. Run the assistant again and answer yes; the summary after the questions lists what was removed. |
| `[x]` "could not generate 3D coordinates" | RDKit could not build a 3D shape for this ligand, usually because of a strained ring system or impossible stereochemistry. Check the structure in your chemistry program, or let the assistant leave it out. |
| Result `CHECK POSE` | The pose did not satisfy the anchoring. See [9.0](#90-reading-the-results-table). |
| Result `FAILED` for a ligand | Read the reason in the table and the `error.txt` in that ligand's folder. Fix the cause, then run `pixi run dock` again with the same output folder and answer yes to resume. |
| `Cannot resume: these inputs or settings differ` | Resuming only works when the files and docking settings are identical to the earlier run. Use a new output folder, or delete the old one. |
| `[x]` "no shared core with the template" | The ligand has nothing in common with the template, so it cannot be anchored. Check it belongs to the same series, or let the assistant leave it out. |
| `[!]` ligand "shares only N heavy atoms" | The ligand is very different from the template. Anchors on such a small core are probably meaningless. |
| `Not enough anchor points found.` | The shared core of a ligand has fewer usable anchor atoms than you asked for. Lower the number of anchors. The assistant normally prevents this by suggesting a safe number. |
| `Failed to load molecule files`, or an error mentioning `NoneType` | RDKit could not read a molecule in one of the SDF files. Open the file in a chemistry program, find the problem molecule and fix or remove it. |
| `File not found` | A typo in the path, or you are not in the `anchored_docking_workflow` folder. Check with `pwd` and `ls`. |
| `You cannot submit to the partition ... : you are not allowed to use the reservation ...` | Your account cannot use that partition or reservation (the H100 nodes need the reservation `h100`, which only some users are on). Choose one of the partitions the assistant lists, usually `ga100`, or ask the cluster administrators for access. The AutoDock-GPU build for large cores works on `ga100` too. |
| `Could not check whether you may use the partition` | The cluster scheduler did not answer. The assistant carries on; if the job is refused later, the reason is printed. |
| Job stays `PD` for a long time | The GPU partition is busy. Wait, or ask for less time (shorter jobs are often scheduled sooner). |
| Job disappeared and there are only some output folders | The time limit ran out. Increase the time limit and rerun (or put the remaining ligands in a new SDF file). |
| Docked ligand is far from the template | The template and receptor are not in the same coordinate system. Re-export both from the **same** structure (see [5.2](#52-the-template-ligand-templatesdf)). |
| Docked ligand clashes with the protein | The receptor may lack hydrogens or contain leftover molecules. Re-check [5.1](#51-the-receptor-receptorpdb). A larger search box or a different seed may also help. |
| `Chain gap between ...` warning | Residues are missing in the structure (an unresolved loop). If the gap is near the template the pocket may be wrong there; the assistant cannot rebuild the loop. Consider a structure that has it, or model it first and use route B. |
| `Protonation with PDB2PQR did not work` | The assistant falls back to ADFR's quicker hydrogens. The message gives the reason. Typical causes are a damaged residue or an unusual kept molecule; try again with that molecule removed. |
| `PDB2PQR: ... hydrogens could not be placed without a clash` or `protonation failed` | That residue, usually next to a gap or a poorly resolved side chain, was left with the best guess PDB2PQR could make. Check it in the protonation report and in PyMOL if it is near the template. |
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
