import os
import shutil
import json
import copy
import argparse
import hashlib
import sys
import traceback
import numpy as np
import subprocess
from typing import NamedTuple
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS

REPO_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO_DIR)
import pose_check  # noqa: E402
import results  # noqa: E402
DEFAULT_ADFR_PATH = os.path.join(REPO_DIR, ".pixi", "envs", "adfr")
ADGPU_BINARY_NAMES = ("adgpu", "autodock_gpu_128wi", "autodock_gpu_64wi", "autodock_gpu_256wi")


def find_autodock_gpu(adg_path=None):
    """Locate the AutoDock-GPU binary.

    If adg_path is given, look there; otherwise search PATH (e.g. after
    `module load GCC/13.2.0 AutoDock-GPU/1.5.3-CUDA`).
    """
    for name in ADGPU_BINARY_NAMES:
        if adg_path:
            candidate = os.path.join(adg_path, name)
            if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                return candidate
        else:
            found = shutil.which(name)
            if found:
                return found
    where = adg_path if adg_path else "PATH"
    raise FileNotFoundError(
        f"AutoDock-GPU binary not found in {where} (tried: {', '.join(ADGPU_BINARY_NAMES)}). "
        "On the HPC run: module load GCC/13.2.0 AutoDock-GPU/1.5.3-CUDA"
    )


def find_autogrid(adg_path=None):
    """Locate autogrid4: the pixi environment first, then adg_path as a fallback."""
    found = shutil.which("autogrid4")
    if found:
        return found
    if adg_path:
        candidate = os.path.join(adg_path, "autogrid4")
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    raise FileNotFoundError("autogrid4 not found. Run this script through `pixi run` so the pixi environment is on PATH.")


class LigandError(Exception):
    """A problem with one ligand. It fails that ligand only; the rest of the batch continues."""


class AnchorSearch(NamedTuple):
    anchor_patterns: list      # (SMARTS, atom index in pattern_mol) for each chosen anchor
    anchor_positions: list     # target coordinates of those anchors, taken from the template
    pattern_mol: object        # the shared core as a molecule
    ptn2template: dict         # core atom index -> template atom index
    pattern_smarts: str        # SMARTS of the shared core


def base36encode(number, length=2):
    chars = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    result = ''
    while number > 0:
        number, remainder = divmod(number, 36)
        result = chars[remainder] + result
    return result.zfill(length)

def calculate_center_of_mass(mol):
    positions = mol.GetConformer().GetPositions()
    masses = [atom.GetMass() for atom in mol.GetAtoms()]
    total_mass = sum(masses)
    weighted_positions = positions * np.array(masses)[:, np.newaxis]
    center_of_mass = np.sum(weighted_positions, axis=0) / total_mass
    return center_of_mass

def get_atype(atom, atom_types):
    mol = atom.GetOwningMol()
    for type_info in atom_types:
        pattern = Chem.MolFromSmarts(type_info['smarts'])
        if pattern is not None:
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                if atom.GetIdx() in match:
                    return type_info['atype']
    return None

def find_equivalent_atoms(mol, atom1_idx, atom2_idx):
    Chem.SanitizeMol(mol, Chem.SANITIZE_PROPERTIES)
    symmetry_classes = Chem.CanonicalRankAtoms(mol, breakTies=False)
    return symmetry_classes[atom1_idx] == symmetry_classes[atom2_idx]

def find_anchors(smart_pattern):
    pattern_mol = Chem.MolFromSmarts(smart_pattern)
    heavy_atoms = [atom.GetIdx() for atom in pattern_mol.GetAtoms() if atom.GetAtomicNum() > 1]
    anchors = []
    for idx in heavy_atoms:
        if all(not find_equivalent_atoms(pattern_mol, idx, idx2) for idx2 in heavy_atoms if idx != idx2):
            anchors.append((smart_pattern, idx))
    return anchors

def spread_selection(positions, count):
    """Indices of up to `count` positions that are as far apart from each other as possible.

    Farthest-point sampling: start with the position farthest from the centre of all of them, then keep
    adding the position whose nearest already-chosen neighbour is farthest away. The result does not
    depend on a random seed. If there are fewer positions than `count`, all of them are returned.
    """
    points = np.asarray(positions, dtype=float)
    count = min(count, len(points))
    if count <= 0:
        return []
    chosen = [int(np.argmax(np.linalg.norm(points - points.mean(axis=0), axis=1)))]
    nearest = np.linalg.norm(points - points[chosen[0]], axis=1)
    while len(chosen) < count:
        nearest[chosen] = -1.0  # never pick the same position twice, even if some coincide
        following = int(np.argmax(nearest))
        chosen.append(following)
        nearest = np.minimum(nearest, np.linalg.norm(points - points[following], axis=1))
    return sorted(chosen)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--template_file", help="Input SDF file", required=True)
    parser.add_argument("--input_file", help="Input SDF file", required=True)
    parser.add_argument("--receptor_file", help="Receptor PDB file", required=True)
    parser.add_argument("--docking_dir", help="Directory for docking", required=True)
    parser.add_argument("--query_smarts", help="Manual SMARTS query for anchor points", default=None)
    parser.add_argument("--anchor_indices", nargs="+", type=int, help="Specify anchor points by atom indices. Query SMILES must be provided. Indices is 0-based.", default=None)
    parser.add_argument("--anchor_num", type=int, default=3, help="Number of anchor points (a maximum in spread mode; 0 means every usable atom)")
    parser.add_argument("--anchor_mode", choices=["random", "area", "spread"], default="random", help="Anchor point selection mode. The number of anchors of area mode is always 3. Spread mode takes up to --anchor_num atoms of the shared core that are as far apart as possible (fewer if the core has fewer usable atoms).")
    parser.add_argument("--size", type=float, help="Size of docking box", default=20.0)
    parser.add_argument("--force", action="store_true", help="Force to overwrite the directory")
    parser.add_argument("--adg_path", default=None,
                        help="Directory containing the AutoDock-GPU binary. Default: search PATH (use `module load GCC/13.2.0 AutoDock-GPU/1.5.3-CUDA`).")
    parser.add_argument("--adfr_path", default=DEFAULT_ADFR_PATH,
                        help="ADFRsuite installation directory. Default: the 'adfr' pixi environment.")
    parser.add_argument("--scripts_path", default="scripts")
    parser.add_argument("--samples", type=int, default=100, help="Number of samples")
    parser.add_argument("--random_seed", type=int, default=None, help="Random seed for sampling")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    parser.add_argument("--clean", action="store_true", help="Clean the directory after docking")
    parser.add_argument("--template_view_file", default=None,
                        help="Only for display: the template with its hydrogens, when --template_file is a heavy-atom "
                             "template whose bond orders are unknown. Without it, missing template hydrogens are added "
                             "with RDKit, which needs real bond orders. Not used for docking.")
    parser.add_argument("--add_hydrogens", action="store_true",
                        help="Let ADFR add hydrogens to the receptor if it has none (a quick guess: histidine "
                             "protonation is not optimised). Receptors that already have hydrogens are untouched.")
    parser.add_argument("--resume", action="store_true",
                        help="Continue a previous run in --docking_dir: ligands that already have a docked pose are "
                             "skipped, the others are redone. Inputs and docking settings must be unchanged.")
    args = parser.parse_args()

    if args.resume and args.force:
        parser.error("--resume keeps the existing results and --force deletes them; use only one.")
    if args.anchor_indices is not None and args.query_smarts is None:
        raise ValueError("Query SMILES must be provided when anchor indices are specified.")
    if args.anchor_indices is not None and args.anchor_num > 0:
        print("Anchor indices are specified. Anchor number will be ignored.")
    if args.anchor_indices is not None and args.anchor_mode in ("area", "spread"):
        print("Anchor indices are specified. Anchor mode will be ignored.")
    return args

SETTINGS_FILE = "run_settings.json"
SUMMARY_FILE = "summary.tsv"
SUMMARY_COLUMNS = ("index", "name", "status", "anchor_max_dev", "core_rmsd", "note")


def file_digest(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def run_signature(template_file, input_file, receptor_file, anchor_num, anchor_mode, size, samples,
                  random_seed, query_smarts, anchor_indices, add_hydrogens=False):
    """Everything that decides what a docking run produces, so a resume can refuse to mix runs."""
    return {
        "add_hydrogens": bool(add_hydrogens),
        "template_sha256": file_digest(template_file),
        "ligands_sha256": file_digest(input_file),
        "receptor_sha256": file_digest(receptor_file),
        "anchor_num": anchor_num,
        "anchor_mode": anchor_mode,
        "size": float(size),
        "samples": samples,
        "random_seed": random_seed,
        "query_smarts": query_smarts,
        "anchor_indices": list(anchor_indices) if anchor_indices else None,
    }


def signature_differences(old, new):
    """Names of the settings that differ between a previous run and the requested one."""
    return [key for key in sorted(set(old) | set(new)) if old.get(key) != new.get(key)]


def read_run_settings(docking_dir):
    """Settings of a previous run, or None if the folder has none."""
    path = os.path.join(docking_dir, SETTINGS_FILE)
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        settings = json.load(fh)
    settings.setdefault("add_hydrogens", False)  # runs saved before this option existed did not add any
    return settings


def write_summary(docking_dir, rows):
    """Write summary.tsv (one row per ligand, in order). Replaced atomically so it is never half-written."""
    path = os.path.join(docking_dir, SUMMARY_FILE)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\t".join(SUMMARY_COLUMNS) + "\n")
        for index in sorted(rows):
            row = rows[index]
            fh.write("\t".join(str(row.get(c, "")).replace("\t", " ").replace("\n", " ") for c in SUMMARY_COLUMNS) + "\n")
    os.replace(tmp, path)


def read_summary(docking_dir):
    """Rows of summary.tsv as a list of dicts (empty if there is none)."""
    path = os.path.join(docking_dir, SUMMARY_FILE)
    if not os.path.isfile(path):
        return []
    with open(path) as fh:
        lines = fh.read().splitlines()
    header = lines[0].split("\t") if lines else []
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]


def ligand_label(mol, index):
    if mol is not None and mol.HasProp("_Name") and mol.GetProp("_Name").strip():
        return mol.GetProp("_Name").strip()
    return f"molecule {index + 1}"


def prepare_docking_directory(docking_dir, force, resume=False, signature=None):
    if resume:
        if not os.path.isdir(docking_dir):
            raise ValueError(f"Nothing to resume: '{docking_dir}' does not exist.")
        previous = read_run_settings(docking_dir)
        if previous is None:
            raise ValueError(f"'{docking_dir}' has no {SETTINGS_FILE}, so it is not a previous run of this "
                             "workflow. Use a new folder, or --force to overwrite.")
        if signature is not None:
            changed = signature_differences(previous, signature)
            if changed:
                raise ValueError("Cannot resume: these inputs or settings differ from the previous run: "
                                 + ", ".join(changed) + ". Use a new folder, or --force to start again.")
        return
    if force:
        if os.path.exists(docking_dir):
            shutil.rmtree(docking_dir)
        os.makedirs(docking_dir)
    else:
        if not os.path.exists(docking_dir):
            os.makedirs(docking_dir)
        else:
            raise ValueError(f"Directory '{docking_dir}' already exists. Use --force to overwrite.")

def prepare_boxsize_file(docking_dir, size, center_of_mass):
    grid_dict = {
        "size_x": size,
        "size_y": size,
        "size_z": size,
        "center_x": center_of_mass[0],
        "center_y": center_of_mass[1],
        "center_z": center_of_mass[2],
    }
    grid_text = "\n".join(f"{key} = {value}" for key, value in grid_dict.items())
    with open(os.path.join(docking_dir, "boxsize.txt"), "w") as f:
        f.write(grid_text)

def find_anchor_points(template_mol, input_mol, anchor_num, samples, random_seed, query_smarts=None, anchor_indices=None, anchor_mode="random", verbose=False):
    if verbose:
        print("Template_SMILES:", Chem.MolToSmiles(template_mol))
        print("Input_SMILES:", Chem.MolToSmiles(input_mol))
    if query_smarts is not None:
        pattern_smiles = query_smarts
        pattern_smarts = query_smarts
        pattern_mol = Chem.MolFromSmiles(pattern_smiles)
        smart_pattern_mol = Chem.MolFromSmarts(pattern_smiles)
        ptn2template_mapping = {
            i: j
            for i, j in zip(
                pattern_mol.GetSubstructMatch(smart_pattern_mol),
                template_mol.GetSubstructMatch(smart_pattern_mol),
            )
        }
    else:
        mcs_res = rdFMCS.FindMCS(
            [input_mol, template_mol],
            atomCompare=rdFMCS.AtomCompare.CompareAny,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            ringCompare=rdFMCS.RingCompare.StrictRingFusion,
            matchValences=False,
            ringMatchesRingOnly=True,
            completeRingsOnly=True,
            matchChiralTag=False,
        )
        substructure = input_mol.GetSubstructMatch(mcs_res.queryMol)
        pattern_smiles = Chem.MolFragmentToSmiles(input_mol, substructure)
        pattern_smarts = Chem.MolToSmarts(Chem.MolFromSmiles(pattern_smiles))
        pattern_mol = Chem.MolFromSmiles(pattern_smiles)
        smart_pattern_mol = Chem.MolFromSmarts(pattern_smarts)
        input2template_mapping = {
            i: j
            for i, j in zip(
                input_mol.GetSubstructMatch(mcs_res.queryMol),
                template_mol.GetSubstructMatch(mcs_res.queryMol),
            )
        }
        ptn2input_mapping = {
            i: j
            for i, j in zip(
                pattern_mol.GetSubstructMatch(smart_pattern_mol),
                input_mol.GetSubstructMatch(smart_pattern_mol),
            )
        }
        ptn2template_mapping = {i: input2template_mapping[j] for i, j in ptn2input_mapping.items()}
    if verbose:
        print("Pattern SMILES:", pattern_smiles)
        print("Pattern SMARTS:", pattern_smarts)

    anchor_patterns_all = find_anchors(pattern_smarts)
    anchor_positions_all = [
        list(template_mol.GetConformer().GetAtomPosition(ptn2template_mapping[idx]))
        for _, idx in anchor_patterns_all
    ]

    if anchor_indices is not None:
        avail_anchors = [idx for anchor,idx in anchor_patterns_all]
        if not all([idx in avail_anchors for idx in anchor_indices]):
            raise ValueError("Some anchor indices are not found in the anchor points, Available anchor indices are: ", avail_anchors)
        anchor_patterns = [(pattern_smarts, idx) for idx in anchor_indices]
        anchor_positions = [
            list(template_mol.GetConformer().GetAtomPosition(ptn2template_mapping[idx]))
            for idx in anchor_indices
        ]
    else:
        if anchor_num > 0:
            rng = np.random.default_rng(random_seed)
            if anchor_mode == "random":
                if len(anchor_patterns_all) < anchor_num:
                    raise ValueError("Not enough anchor points found.")
                sample_idx = rng.choice(len(anchor_patterns_all), size=min(anchor_num, len(anchor_patterns_all)), replace=False)
                anchor_patterns = [anchor_patterns_all[idx] for idx in sample_idx]
                anchor_positions = [anchor_positions_all[idx] for idx in sample_idx]
            elif anchor_mode == "spread":
                chosen = spread_selection(anchor_positions_all, anchor_num)
                anchor_patterns = [anchor_patterns_all[idx] for idx in chosen]
                anchor_positions = [anchor_positions_all[idx] for idx in chosen]
            elif anchor_mode == "area":
                anchor_num = 3
                area_max = -1
                anchor_patterns = None
                anchor_positions = None
                for _ in range(samples):
                    if len(anchor_patterns_all) < anchor_num:
                        raise ValueError("Not enough anchor points found.")
                    sample_idx = rng.choice(len(anchor_patterns_all), size=anchor_num, replace=False)
                    positions = [np.array(anchor_positions_all[idx]) for idx in sample_idx]
                    if len(positions) != 3:
                        raise ValueError("Not enough anchor points found.")
                    area = 0.5 * np.linalg.norm(
                        np.cross(positions[1] - positions[0], positions[2] - positions[0])
                    )
                    if area > area_max:
                        area_max = area
                        anchor_patterns = [anchor_patterns_all[idx] for idx in sample_idx]
                        anchor_positions = [anchor_positions_all[idx] for idx in sample_idx]
            else:
                raise ValueError(f"Unknown anchor mode: {anchor_mode}")
        else:
            anchor_patterns = anchor_patterns_all
            anchor_positions = anchor_positions_all
    return AnchorSearch(anchor_patterns, anchor_positions, pattern_mol, ptn2template_mapping, pattern_smarts)


def generate_parameters_json(scripts_path, anchor_patterns, pattern_mol):
    with open(os.path.join(scripts_path, "base_parameters.json")) as f:
        base_parameters = json.load(f)
    atom_types = base_parameters['ad4']
    parameters = copy.deepcopy(base_parameters)

    new_atype2atype = {}
    anchor_dict_list = []

    for i, (smarts, idx) in enumerate(anchor_patterns):
        atom = pattern_mol.GetAtomWithIdx(idx)
        atype = get_atype(atom, atom_types)
        if atype is None:
            raise ValueError(f"Could not determine atom type for atom index {idx}.")
        new_atype = 'A' + base36encode(i, 2)
        d = {"smarts": smarts, "IDX": [idx], "atype": new_atype}
        parameters['ad4'].append(d)
        anchor_dict_list.append(d)
        new_atype2atype[new_atype] = atype
    return parameters, anchor_dict_list, new_atype2atype

# Each embedding attempt is cut off after this long. Ordinary ligands take well under a second;
# a hopeless one (measured: a strained cage) otherwise needs about a minute to give up.
EMBED_TIMEOUT_SECONDS = 30


def embed_ligand(input_mol, timeout=EMBED_TIMEOUT_SECONDS):
    """Add hydrogens and generate 3D coordinates.

    RDKit signals a failed embedding only through its return value, and leaves the molecule
    without any coordinates; writing it out then produces a flat 2D layout that Meeko accepts.
    Raise instead, so such a ligand is reported and not docked with a wrong geometry.
    """
    mol_h = Chem.AddHs(input_mol, addCoords=True)
    params = AllChem.ETKDGv3()
    params.timeout = timeout
    if AllChem.EmbedMolecule(mol_h, params) != 0:
        retry = AllChem.ETKDGv3()
        retry.timeout = timeout
        retry.useRandomCoords = True  # RDKit's usual fallback for hard-to-embed molecules
        if AllChem.EmbedMolecule(mol_h, retry) != 0:
            raise LigandError("RDKit could not generate 3D coordinates for this ligand "
                              "(strained or impossible geometry or stereochemistry)")
    return mol_h

def prepare_ligand(docking_dir, input_mol):
    inputmol_h = embed_ligand(input_mol)
    ligand_sdf = os.path.join(docking_dir, "ligand.sdf")
    with Chem.SDWriter(ligand_sdf) as f:
        f.write(inputmol_h)
    subprocess.run(
        ["mk_prepare_ligand.py", "-p", "parameters.json", "-i", "ligand.sdf", "-o", "ligand.pdbqt"],
        cwd=docking_dir,
        check=True
    )

def prepare_receptor(docking_dir, receptor_file, adfr_path, add_hydrogens=False):
    receptor_pdb = os.path.join(docking_dir, "receptor.pdb")
    shutil.copy(receptor_file, receptor_pdb)
    prepare_receptor_cmd = [
        os.path.join(adfr_path, "bin", "prepare_receptor"),
        "-r",
        "receptor.pdb",
        "-o",
        "receptor.pdbqt",
    ]
    if add_hydrogens:
        prepare_receptor_cmd += ["-A", "checkhydrogens"]  # only adds hydrogens if the file has none
    subprocess.run(prepare_receptor_cmd, cwd=docking_dir, check=True)

def prepare_grid(docking_dir, scripts_path, autogrid_bin):
    write_gpf_cmd = [
        os.path.join(scripts_path, "write-gpf.py"),
        "receptor.pdbqt",
        "-b",
        "boxsize.txt",
        "--mapprefix",
        "rec",
    ]
    subprocess.run(write_gpf_cmd, cwd=docking_dir, check=True)
    autogrid_cmd = [autogrid_bin, "-p", "rec.gpf", "-l", "rec.glg"]
    subprocess.run(autogrid_cmd, cwd=docking_dir, check=True)

def add_bias_to_map_file(docking_dir, scripts_path, anchor_dict_list, anchor_positions, new_atype2atype):
    for anchor_d, position in zip(anchor_dict_list, anchor_positions):
        atype = anchor_d["atype"]
        old_atype = new_atype2atype[atype]
        addbias_cmd = [
            os.path.join(scripts_path, "addbias.py"),
            "-i",
            f"rec.{old_atype}.map",
            "-o",
            f"rec.{atype}.map",
            "-x",
            str(position[0]),
            str(position[1]),
            str(position[2]),
        ]
        subprocess.run(addbias_cmd, cwd=docking_dir, check=True)
        insert_type_cmd = [
            os.path.join(scripts_path, "insert_type_in_fld.py"),
            "rec.maps.fld",
            "--newtype",
            atype,
        ]
        subprocess.run(insert_type_cmd, cwd=docking_dir, check=True)

def run_docking(docking_dir, adgpu_bin, anchor_dict_list, new_atype2atype):
    t_options_map = {}
    for anchor_d in anchor_dict_list:
        atype = anchor_d["atype"]
        old_atype = new_atype2atype[atype]
        t_options_map.setdefault(old_atype, []).append(atype)
    t_options = "/".join([f"{','.join(v)}={k}" for k, v in t_options_map.items()])

    adgpu_cmd = [
        adgpu_bin,
        "-L",
        "ligand.pdbqt",
        "-M",
        "rec.maps.fld",
        "-T",
        t_options,
    ]
    subprocess.run(adgpu_cmd, cwd=docking_dir, check=True)

def export_results(docking_dir):
    export_cmd = ["mk_export.py", "ligand.dlg", "-s", "ligand_docked.sdf"]
    subprocess.run(export_cmd, cwd=docking_dir, check=True)

def generate_pymol_session(docking_dir, template_file, receptor_file):
    """Write a PyMOL script that loads the receptor, template, and docked ligand."""
    template_sdf = os.path.join(docking_dir, "template.sdf")
    shutil.copy2(template_file, template_sdf)

    pml_content = """\
# Load structures. Run this script from inside this folder (pymol view_results.pml).
# The receptor and template are the versions with every hydrogen, one folder up.
load ../receptor_with_hydrogens.pdb, receptor
load ../template_with_hydrogens.sdf, template
load ligand_docked.sdf, docked_ligand

# Style receptor
hide everything, receptor
show cartoon, receptor
color gray80, receptor
set cartoon_transparency, 0.5, receptor

# Style template ligand
show sticks, template
color cyan, template

# Style docked ligand
show sticks, docked_ligand
color magenta, docked_ligand

# General settings
bg_color white
set stick_radius, 0.15
set ray_opaque_background, 1
zoom docked_ligand, 5

# Save session
save docking_results.pse
"""
    pml_path = os.path.join(docking_dir, "view_results.pml")
    with open(pml_path, "w") as f:
        f.write(pml_content)


def clean_directory(docking_dir):
    keep = {"ligand_docked.sdf", "template.sdf", "receptor.pdb", "view_results.pml", "docking_results.pse",
            "pose_check.json", "anchor_positions.json"}
    for f in os.listdir(docking_dir):
        if f not in keep:
            os.remove(os.path.join(docking_dir, f))

def main():
    args = parse_arguments()

    adg_path = os.path.abspath(args.adg_path) if args.adg_path else None
    adfr_path = os.path.abspath(args.adfr_path)
    scripts_path = os.path.abspath(args.scripts_path)
    adgpu_bin = find_autodock_gpu(adg_path)
    autogrid_bin = find_autogrid(adg_path)

    signature = run_signature(args.template_file, args.input_file, args.receptor_file, args.anchor_num,
                              args.anchor_mode, args.size, args.samples, args.random_seed,
                              args.query_smarts, args.anchor_indices, args.add_hydrogens)

    # Prepare docking directory
    prepare_docking_directory(args.docking_dir, args.force, args.resume, signature)
    if not args.resume:
        with open(os.path.join(args.docking_dir, SETTINGS_FILE), "w") as f:
            json.dump(signature, f, indent=4)

    # Load molecules
    template_mol = Chem.SDMolSupplier(args.template_file)[0]
    input_mols = list(Chem.SDMolSupplier(args.input_file))
    if template_mol is None or input_mols is None:
        raise ValueError("Failed to load molecule files.")

    tools = {"adfr_path": adfr_path, "scripts_path": scripts_path, "adgpu_bin": adgpu_bin,
             "autogrid_bin": autogrid_bin}
    prepare_view_files(args, adfr_path)
    rows = {}
    total = len(input_mols)
    for index, input_mol in enumerate(input_mols):
        name = ligand_label(input_mol, index)
        docking_dir = os.path.join(args.docking_dir, str(index))
        pose_path = os.path.join(docking_dir, "ligand_docked.sdf")

        if args.resume and os.path.isfile(pose_path):
            rows[index] = previous_row(index, name, docking_dir)
            print(f"[{index + 1}/{total}] {name}: already docked, skipped", flush=True)
            write_summary(args.docking_dir, rows)
            continue

        if os.path.isdir(docking_dir):  # left over from an interrupted or failed attempt
            shutil.rmtree(docking_dir)
        os.makedirs(docking_dir)
        try:
            if input_mol is None:
                raise LigandError("RDKit could not read this molecule from the input file")
            checked = dock_ligand(args, tools, template_mol, input_mol, docking_dir)
        except Exception as exc:  # one ligand must not stop the others
            note = describe_failure(exc)
            with open(os.path.join(docking_dir, "error.txt"), "w") as f:
                f.write(f"{note}\n\n{traceback.format_exc()}")
            rows[index] = {"index": index, "name": name, "status": "failed", "note": note}
            print(f"[{index + 1}/{total}] {name}: FAILED - {note}", flush=True)
        else:
            rows[index] = row_from_check(index, name, checked)
            print(f"[{index + 1}/{total}] {name}: {describe_row(rows[index])}", flush=True)
            if args.clean:
                clean_directory(docking_dir)
        write_summary(args.docking_dir, rows)

    try:
        made = results.write_all(args.docking_dir)
        print(f"Combined results: {results.RESULTS_CSV}, {results.TOP_POSES}, {results.ALL_POSES} "
              f"({made['poses']} poses) and {results.VIEW_SCRIPT} in {args.docking_dir}")
    except Exception as exc:  # the docking is done; a problem with the summary files must not hide that
        print(f"Could not write the combined result files ({type(exc).__name__}: {exc}). "
              f"You can retry with: python results.py {args.docking_dir}")

    failed = [r for r in rows.values() if r["status"] == "failed"]
    flagged = [r for r in rows.values() if r["status"] == "flagged"]
    docked = [r for r in rows.values() if r["status"] in ("docked", "flagged")]
    print(f"\nFinished: {len(docked)} of {total} ligand(s) docked, {len(failed)} failed, "
          f"{len(flagged)} pose(s) flagged for checking.")
    if flagged:
        print("Flagged poses did not satisfy the anchoring; look at them before trusting them "
              f"(see {SUMMARY_FILE}).")
    for r in failed:
        print(f"  ligand {r['index']} ({r['name']}): {r['note']}")
    if failed:
        print(f"Details are in <docking_dir>/<index>/error.txt. Fix the cause and run again with --resume "
              "to redo only the failed ligands.")
    print("To visualise a result, open its PyMOL script: pymol <docking_dir>/<index>/view_results.pml")
    if failed:
        sys.exit(1)


def prepare_view_files(args, adfr_path):
    """Write the receptor and template with every hydrogen, for looking at the results.

    AutoDock's own files have the non-polar hydrogens merged away, so they are not what to look at.
    A failure here is reported but never stops the docking.
    """
    out = args.docking_dir
    try:
        n, how = results.make_view_receptor(args.receptor_file, adfr_path, args.add_hydrogens,
                                            os.path.join(out, results.VIEW_RECEPTOR))
        print(f"Receptor for viewing: {n} hydrogens ({how})", flush=True)
    except Exception as exc:
        print(f"Could not prepare the receptor with hydrogens for viewing ({type(exc).__name__}: {exc}); "
              "the unmodified receptor is used.", flush=True)
        shutil.copyfile(args.receptor_file, os.path.join(out, results.VIEW_RECEPTOR))
    try:
        n, how = results.make_view_template(args.template_file, os.path.join(out, results.VIEW_TEMPLATE),
                                            args.template_view_file)
        print(f"Template for viewing: {n} hydrogens ({how})", flush=True)
    except Exception as exc:
        print(f"Could not prepare the template with hydrogens for viewing ({type(exc).__name__}: {exc}); "
              "the unmodified template is used.", flush=True)
        shutil.copyfile(args.template_file, os.path.join(out, results.VIEW_TEMPLATE))


def dock_ligand(args, tools, template_mol, input_mol, docking_dir):
    """Run the whole workflow for one ligand. Returns the pose check result."""
    # Calculate center of mass for the grid center
    center_of_mass = calculate_center_of_mass(template_mol)

    # Prepare boxsize.txt
    prepare_boxsize_file(docking_dir, args.size, center_of_mass)

    # Find anchor points
    search = find_anchor_points(
        template_mol,
        input_mol,
        args.anchor_num,
        args.samples,
        args.random_seed,
        args.query_smarts,
        args.anchor_indices,
        args.anchor_mode,
        args.verbose,
    )
    with open(os.path.join(docking_dir, "anchor_positions.json"), "w") as f:
        json.dump(search.anchor_positions, f, indent=4)

    # Generate parameters.json
    parameters, anchor_dict_list, new_atype2atype = generate_parameters_json(
        tools["scripts_path"],
        search.anchor_patterns,
        search.pattern_mol,
    )

    # Save parameters.json
    with open(os.path.join(docking_dir, "parameters.json"), "w") as f:
        json.dump(parameters, f, indent=4)

    # Prepare ligand.pdbqt
    prepare_ligand(docking_dir, input_mol)

    # Prepare receptor.pdbqt
    prepare_receptor(docking_dir, args.receptor_file, tools["adfr_path"], args.add_hydrogens)

    # Prepare grid
    prepare_grid(docking_dir, tools["scripts_path"], tools["autogrid_bin"])

    # Add bias to the map file
    add_bias_to_map_file(
        docking_dir,
        tools["scripts_path"],
        anchor_dict_list,
        search.anchor_positions,
        new_atype2atype,
    )

    # Run docking
    run_docking(docking_dir, tools["adgpu_bin"], anchor_dict_list, new_atype2atype)

    # Export results
    export_results(docking_dir)

    # Generate PyMOL session script
    generate_pymol_session(docking_dir, args.template_file, args.receptor_file)

    # Check the pose against the anchoring
    return check_docked_pose(docking_dir, template_mol, search)


def check_docked_pose(docking_dir, template_mol, search):
    """Measure how well the top pose satisfies the anchoring and save it as pose_check.json."""
    try:
        checked = pose_check.check_docked_file(
            os.path.join(docking_dir, "ligand_docked.sdf"), template_mol, search.pattern_mol,
            search.pattern_smarts, search.ptn2template, [idx for _, idx in search.anchor_patterns],
            search.anchor_positions)
    except Exception as exc:  # a failed check must not discard a finished docking
        checked = {"ok": None, "problems": [f"the pose check could not be run ({type(exc).__name__}: {exc})"],
                   "anchor_max": None, "core_rmsd": None, "anchor_deviations": [], "core_atoms": 0}
    checked["anchor_tolerance"] = pose_check.ANCHOR_TOLERANCE
    checked["core_rmsd_tolerance"] = pose_check.CORE_RMSD_TOLERANCE
    with open(os.path.join(docking_dir, "pose_check.json"), "w") as f:
        json.dump(checked, f, indent=4)
    return checked


def row_from_check(index, name, checked):
    status = "flagged" if checked.get("ok") is False else "docked"
    return {"index": index, "name": name, "status": status,
            "anchor_max_dev": "" if checked.get("anchor_max") is None else checked["anchor_max"],
            "core_rmsd": "" if checked.get("core_rmsd") is None else checked["core_rmsd"],
            "note": "; ".join(checked.get("problems", []))}


def previous_row(index, name, docking_dir):
    """Summary row for a ligand docked in an earlier run, from the pose check saved then."""
    path = os.path.join(docking_dir, "pose_check.json")
    if os.path.isfile(path):
        with open(path) as f:
            return row_from_check(index, name, json.load(f))
    return {"index": index, "name": name, "status": "docked", "note": "pose was not checked"}


def describe_row(row):
    if row["status"] == "flagged":
        return f"docked, but flagged: {row['note']}"
    text = "docked"
    if row.get("anchor_max_dev") not in ("", None):
        text += f" (anchor deviation {row['anchor_max_dev']} A, core RMSD {row.get('core_rmsd', '')} A)"
    elif row.get("note"):
        text += f" ({row['note']})"
    return text


def describe_failure(exc):
    """One line saying what went wrong for a ligand."""
    if isinstance(exc, subprocess.CalledProcessError):
        program = os.path.basename(exc.cmd[0]) if isinstance(exc.cmd, (list, tuple)) and exc.cmd else str(exc.cmd)
        return f"{program} failed (exit code {exc.returncode})"
    text = str(exc).strip().splitlines()[0] if str(exc).strip() else type(exc).__name__
    return text if isinstance(exc, LigandError) else f"{type(exc).__name__}: {text}"[:300]


if __name__ == "__main__":
    main()
