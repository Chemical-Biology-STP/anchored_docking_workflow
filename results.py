"""Combined outputs for a docking run, and the files needed to look at it with every hydrogen present.

AutoDock works with non-polar hydrogens merged into their carbon, so its own files are not what you
want to look at. The pieces here:

- the docked ligand: Meeko already writes the complete molecule (checked: every hydrogen is there, with
  sensible X-H lengths), so the poses are used as they are;
- the receptor: ADFR's docking input has only the polar hydrogens, so the display copy is prepared again
  without merging. Its polar hydrogens are at exactly the coordinates used for docking (checked);
- the template: kept as a separate file with hydrogens, for display only.

Run on its own to rebuild the combined files:  pixi run python results.py <docking_dir>
"""

import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

VIEW_RECEPTOR = "receptor_with_hydrogens.pdb"
VIEW_TEMPLATE = "template_with_hydrogens.sdf"
RESULTS_CSV = "results.csv"
TOP_POSES = "top_poses.sdf"
ALL_POSES = "all_poses.sdf"
VIEW_SCRIPT = "view_all.pml"
SESSION_NAME = "docking_all.pse"

CSV_COLUMNS = ("index", "name", "status", "free_energy", "intermolecular_energy", "internal_energy",
               "cluster_size", "n_poses", "anchor_max_dev", "core_rmsd", "folder", "note")

# AutoDock atom types that are not simply the element symbol.
_ELEMENT_OF_TYPE = {"A": "C", "HD": "H", "NA": "N", "OA": "O", "SA": "S"}


# --------------------------------------------------------------------------
# Hydrogens for display
# --------------------------------------------------------------------------


def count_hydrogens(pdb_path):
    """Hydrogen atoms in a PDB file (by element column, or by atom name if there is none)."""
    n = 0
    with open(pdb_path, errors="replace") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                element = line[76:78].strip().upper() or re.sub(r"[^A-Za-z]", "", line[12:16].strip())[:1].upper()
                n += element in ("H", "D")
    return n


def pdbqt_to_pdb(text):
    """PDB text from a receptor PDBQT: coordinates and names kept, AutoDock charges and types dropped.

    PDBQT atom types are not elements (A is aromatic carbon, HD a polar hydrogen, and so on), so the
    element column is rebuilt from the type.
    """
    out = []
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        ad_type = line[77:79].strip() or line.split()[-1]
        element = _ELEMENT_OF_TYPE.get(ad_type, ad_type).upper()
        out.append(f"{line[:54]}  1.00  0.00          {element:>2}")
    out.append("END")
    return "\n".join(out) + "\n"


def make_view_receptor(receptor_file, adfr_path, add_hydrogens, destination):
    """Write the receptor with every hydrogen to `destination`. Returns (hydrogens, how it was made).

    - If the supplied receptor already has hydrogens it is used as it is: ADFR keeps its polar hydrogens
      in the docking input, so this is the receptor that was docked with its non-polar hydrogens put back.
    - If it has none and ADFR added them for docking, ADFR is run again without merging the non-polar
      ones; the polar hydrogens come out at identical coordinates.
    - Otherwise there are no hydrogens to show.
    """
    present = count_hydrogens(receptor_file)
    if present:
        shutil.copyfile(receptor_file, destination)
        return present, "hydrogens as supplied in the receptor file"
    if not add_hydrogens:
        shutil.copyfile(receptor_file, destination)
        return 0, "the receptor has no hydrogens and none were added"
    with tempfile.TemporaryDirectory(prefix="view_receptor_") as folder:
        shutil.copyfile(receptor_file, os.path.join(folder, "receptor.pdb"))
        command = [os.path.join(adfr_path, "bin", "prepare_receptor"), "-r", "receptor.pdb", "-o", "full.pdbqt",
                   "-A", "checkhydrogens", "-U", "deleteAltB"]  # deleteAltB only: keeps the non-polar hydrogens
        subprocess.run(command, cwd=folder, check=True, capture_output=True)
        with open(os.path.join(folder, "full.pdbqt")) as fh:
            text = pdbqt_to_pdb(fh.read())
    with open(destination, "w") as fh:
        fh.write(text)
    return count_hydrogens(destination), "hydrogens added by ADFR, including the non-polar ones"


def make_view_template(template_file, destination, view_file=None):
    """Write the template with hydrogens to `destination`. Returns (hydrogens, how it was made).

    view_file: a ready-made display version (for a template taken from a structure file, whose bond
    orders are unknown). Without it the template is assumed to have real bond orders, and any missing
    hydrogens are added.
    """
    if view_file:
        shutil.copyfile(view_file, destination)
        mol = Chem.SDMolSupplier(destination, removeHs=False)[0]
        n = sum(a.GetAtomicNum() == 1 for a in mol.GetAtoms()) if mol is not None else 0
        return n, ("the prepared display version" if n else "the template has no hydrogens (none could be derived)")
    mol = Chem.SDMolSupplier(template_file, removeHs=False)[0]
    if mol is None:
        shutil.copyfile(template_file, destination)
        return 0, "the template could not be read for display"
    with_h = Chem.AddHs(mol, addCoords=True)  # adds only the missing ones
    with Chem.SDWriter(destination) as writer:
        writer.write(with_h)
    return sum(a.GetAtomicNum() == 1 for a in with_h.GetAtoms()), "hydrogens added with RDKit"


# --------------------------------------------------------------------------
# Scores and combined files
# --------------------------------------------------------------------------


def pose_scores(mol):
    """Scores Meeko stored on a docked pose (a JSON string in the 'meeko' property)."""
    if mol is None or not mol.HasProp("meeko"):
        return {}
    try:
        return json.loads(mol.GetProp("meeko"))
    except ValueError:
        return {}


def read_poses(path):
    try:
        return [m for m in Chem.SDMolSupplier(path, removeHs=False) if m is not None]
    except OSError:
        return []


def top_pose_scores(path):
    """Scores of the best pose in a docked SDF, plus how many poses it holds."""
    poses = read_poses(path)
    scores = pose_scores(poses[0]) if poses else {}
    return {"free_energy": scores.get("free_energy", ""), "intermolecular_energy": scores.get("intermolecular_energy", ""),
            "internal_energy": scores.get("internal_energy", ""), "cluster_size": scores.get("cluster_size", ""),
            "n_poses": len(poses)}


def read_summary(docking_dir):
    path = os.path.join(docking_dir, "summary.tsv")
    if not os.path.isfile(path):
        return []
    with open(path) as fh:
        lines = fh.read().splitlines()
    header = lines[0].split("\t") if lines else []
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]


def collect(docking_dir):
    """One row per ligand: status and pose-check numbers from summary.tsv, scores from the docked pose."""
    rows = []
    for row in read_summary(docking_dir):
        index = row.get("index", "")
        pose = os.path.join(docking_dir, index, "ligand_docked.sdf")
        entry = {"index": index, "name": row.get("name", ""), "status": row.get("status", ""),
                 "anchor_max_dev": row.get("anchor_max_dev", ""), "core_rmsd": row.get("core_rmsd", ""),
                 "folder": index, "note": row.get("note", "")}
        if os.path.isfile(pose):
            entry.update(top_pose_scores(pose))
        else:  # a failed ligand has no pose: keep the same columns, empty
            entry.update({"free_energy": "", "intermolecular_energy": "", "internal_energy": "",
                          "cluster_size": "", "n_poses": 0})
        rows.append(entry)
    return rows


def write_csv(docking_dir, rows):
    path = os.path.join(docking_dir, RESULTS_CSV)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def write_pose_files(docking_dir, rows):
    """top_poses.sdf (best pose of each docked ligand) and all_poses.sdf (every cluster lead).

    Each record is the complete molecule, with hydrogens, and carries the scores as data fields.
    """
    top = Chem.SDWriter(os.path.join(docking_dir, TOP_POSES))
    everything = Chem.SDWriter(os.path.join(docking_dir, ALL_POSES))
    written = 0
    try:
        for row in rows:
            if row["status"] not in ("docked", "flagged"):
                continue
            poses = read_poses(os.path.join(docking_dir, row["folder"], "ligand_docked.sdf"))
            for rank, pose in enumerate(poses, 1):
                scores = pose_scores(pose)
                pose.SetProp("_Name", row["name"])
                pose.SetProp("ligand_index", str(row["index"]))
                pose.SetProp("pose_rank", str(rank))
                for key in ("free_energy", "intermolecular_energy", "internal_energy", "cluster_size"):
                    if key in scores:
                        pose.SetProp(key, str(scores[key]))
                if rank == 1:  # the pose check is made on the best pose only
                    pose.SetProp("status", row["status"])
                    for key in ("anchor_max_dev", "core_rmsd"):
                        if row.get(key) not in ("", None):
                            pose.SetProp(key, str(row[key]))
                    top.write(pose)
                everything.write(pose)
                written += 1
    finally:
        top.close()
        everything.close()
    return written


def pymol_name(index, name):
    return f"L{index}_" + re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")[:30]


def write_view_script(docking_dir, rows):
    """A PyMOL script that loads the protein, the template and every docked ligand into one session."""
    docked = [r for r in rows if r["status"] in ("docked", "flagged")]
    lines = ["# All docked ligands in one PyMOL session. Run it from this folder:",
             f"#   pymol {VIEW_SCRIPT}      (or, inside PyMOL:  @{VIEW_SCRIPT})",
             f"load {VIEW_RECEPTOR}, receptor", f"load {VIEW_TEMPLATE}, template", ""]
    for row in docked:
        lines.append(f"load {row['folder']}/ligand_docked.sdf, {pymol_name(row['index'], row['name'])}")
    names = [pymol_name(r["index"], r["name"]) for r in docked]
    flagged = [pymol_name(r["index"], r["name"]) for r in docked if r["status"] == "flagged"]
    ok_names = [n for n in names if n not in flagged]
    lines += ["", "# Style: protein as cartoon with the binding site shown as sticks, hydrogens included",
              "bg_color white", "hide everything, receptor", "show cartoon, receptor", "color gray80, receptor",
              "set cartoon_transparency, 0.5, receptor", "set stick_radius, 0.15",
              "show sticks, template", "color cyan, template"]
    if names:
        lines += [f"group ligands, {' '.join(names)}", "show sticks, ligands", "color magenta, ligands and elem C"]
        if flagged:
            lines += ["# Orange: the pose check flagged these poses, look at them before trusting them",
                      f"color orange, ({' or '.join(flagged)}) and elem C"]
        lines += ["select site, byres (receptor within 5 of ligands)", "show sticks, site",
                  "color gray60, site and elem C", "zoom ligands, 6"]
    else:
        lines.append("zoom template, 8")
    lines += ["", "# Ligands with several poses hold them as states: use the arrows at the bottom right.",
              f"save {SESSION_NAME}"]
    path = os.path.join(docking_dir, VIEW_SCRIPT)
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


def write_all(docking_dir):
    """Write every combined file for the run in `docking_dir`. Returns a dict describing what was made."""
    rows = collect(docking_dir)
    write_csv(docking_dir, rows)
    poses = write_pose_files(docking_dir, rows)
    write_view_script(docking_dir, rows)
    return {"ligands": len(rows), "docked": sum(r["status"] in ("docked", "flagged") for r in rows), "poses": poses}


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1 or not os.path.isdir(argv[0]):
        print("usage: python results.py <docking_dir>")
        return 2
    info = write_all(argv[0])
    print(f"{info['docked']} of {info['ligands']} ligand(s) docked; wrote {RESULTS_CSV}, {TOP_POSES}, "
          f"{ALL_POSES} ({info['poses']} poses) and {VIEW_SCRIPT} in {argv[0]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
