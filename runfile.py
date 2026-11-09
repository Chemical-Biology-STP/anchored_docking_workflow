"""The run file: a YAML description of one docking run, used as `dock run.yaml` instead of the interactive wizard.

This module only reads, checks and writes the file. It knows nothing about structures or ligands; the checks that
need the data (which molecule is the template, are there flagged ligands, ...) are made by dock_cli.py, which stops
with an explanation when the file does not answer a question the wizard would have asked.
"""
import difflib
import os
import re

import yaml

STRUCTURE_EXTENSIONS = (".pdb", ".ent", ".cif", ".mmcif")

TOP_KEYS = ("structure", "template_molecule", "chains", "keep_molecules", "receptor", "template", "ligands",
            "output", "existing_output", "protonation", "ph", "flagged_ligands", "docking", "run")
DOCKING_KEYS = ("anchor_mode", "max_anchors", "box_size", "samples", "seed", "core_smarts", "anchor_atoms",
                "clean_intermediates", "verbose")
RUN_KEYS = ("mode", "job_name", "time", "partition", "reservation", "gpus", "modules", "pixi_dir",
            "autodock_gpu", "watch")


class RunFileError(Exception):
    """The run file cannot be used. The message says what is wrong and what to do about it."""


# --------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------


class _Problems:
    def __init__(self):
        self.lines = []

    def add(self, key, message):
        self.lines.append(f"{key}: {message}")


def _unknown(keys, allowed, where, problems):
    for key in keys:
        if key not in allowed:
            close = difflib.get_close_matches(str(key), allowed, n=1)
            hint = f" Did you mean '{close[0]}'?" if close else ""
            problems.add(where + str(key), f"unknown key.{hint}")


def _text(value, key, problems):
    if isinstance(value, str) and value.strip():
        return value.strip()
    problems.add(key, "must be text")
    return None


def _path(value, key, base_dir, problems, extensions=(), must_exist=True):
    text = _text(value, key, problems)
    if text is None:
        return None
    path = os.path.abspath(os.path.join(base_dir, os.path.expanduser(os.path.expandvars(text))))
    if extensions and not path.lower().endswith(extensions):
        problems.add(key, f"expected a {' or '.join(extensions)} file, got {os.path.basename(path)}")
    elif must_exist and not os.path.isfile(path):
        problems.add(key, f"file not found: {path}")
    elif must_exist and os.path.getsize(path) == 0:
        problems.add(key, f"the file is empty: {path}")
    return path


def _choice(value, key, choices, problems):
    if isinstance(value, str) and value.lower() in choices:
        return value.lower()
    problems.add(key, f"must be one of: {', '.join(choices)}")
    return None


def _number(value, key, problems, cast=float, minimum=None, maximum=None):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or (cast is int and value != int(value)):
        problems.add(key, "must be a number" if cast is float else "must be a whole number")
        return None
    value = cast(value)
    if minimum is not None and value < minimum:
        problems.add(key, f"must be at least {minimum}")
        return None
    if maximum is not None and value > maximum:
        problems.add(key, f"must be at most {maximum}")
        return None
    return value


def _flag(value, key, problems):
    if isinstance(value, bool):
        return value
    problems.add(key, "must be true or false")
    return None


def _names(value, key, problems):
    """A list of names given as a YAML list or as one text ("A B" or "A, B")."""
    if isinstance(value, str):
        value = value.replace(",", " ").split()
    if isinstance(value, list) and all(isinstance(v, (str, int)) and not isinstance(v, bool) for v in value):
        return [str(v).strip() for v in value if str(v).strip()]
    problems.add(key, "must be a list of names, for example [A, B]")
    return None


def _section(data, name, allowed, problems):
    section = data.get(name)
    if section is None:
        return {}
    if not isinstance(section, dict):
        problems.add(name, "must be a group of settings (key: value lines indented below it)")
        return {}
    _unknown(section, allowed, f"{name}.", problems)
    return section


def load(path):
    """Read a run file and check everything that can be checked without looking at the molecules.

    Returns a dict of the settings given (paths absolute; keys that are not in the file are missing from it).
    Raises RunFileError listing every problem found, not only the first.
    """
    path = os.path.abspath(path)
    base = os.path.dirname(path)
    try:
        with open(path) as fh:
            data = yaml.safe_load(fh)
    except OSError as exc:
        raise RunFileError(f"Cannot read {path}: {exc.strerror or exc}")
    except yaml.YAMLError as exc:
        raise RunFileError(f"{path} is not valid YAML:\n{exc}")
    if data is None:
        data = {}
    if not isinstance(data, dict):
        raise RunFileError(f"{path} must contain 'key: value' lines.")

    problems = _Problems()
    _unknown(data, TOP_KEYS, "", problems)
    cfg = {"file": path, "base_dir": base}

    # Where the protein comes from: a complex to clean up, or a receptor and template that are ready
    has_structure = data.get("structure") is not None
    has_files = data.get("receptor") is not None or data.get("template") is not None
    if has_structure and has_files:
        problems.add("structure", "give either 'structure' (a complex to clean up) or 'receptor' and 'template' "
                                  "(files that are ready to use), not both")
    elif not has_structure and not has_files:
        problems.add("structure", "say what to dock into: either 'structure' (a PDB or mmCIF file of a protein-ligand "
                                  "complex, cleaned up for you) or both 'receptor' (PDB) and 'template' (SDF)")
    elif has_structure:
        cfg["route"] = "structure"
        cfg["structure"] = _path(data["structure"], "structure", base, problems, must_exist=False)
        if cfg["structure"] and not cfg["structure"].lower().endswith(STRUCTURE_EXTENSIONS):
            problems.add("structure", f"expected a {' or '.join(STRUCTURE_EXTENSIONS)} file")
        elif cfg["structure"] and not os.path.isfile(cfg["structure"]):
            problems.add("structure", f"file not found: {cfg['structure']}")
        if data.get("template_molecule") is not None:
            cfg["template_molecule"] = _text(data["template_molecule"], "template_molecule", problems)
        if data.get("chains") is not None:
            cfg["chains"] = _names(data["chains"], "chains", problems)
        if data.get("keep_molecules") is not None:
            cfg["keep_molecules"] = _names(data["keep_molecules"], "keep_molecules", problems)
    else:
        cfg["route"] = "files"
        if data.get("receptor") is None:
            problems.add("receptor", "missing. With a ready 'template' file you also need the 'receptor' PDB file")
        else:
            cfg["receptor"] = _path(data["receptor"], "receptor", base, problems, (".pdb",))
        if data.get("template") is None:
            problems.add("template", "missing. With a 'receptor' file you also need the 'template' SDF file")
        else:
            cfg["template"] = _path(data["template"], "template", base, problems, (".sdf",))
        for key in ("template_molecule", "chains", "keep_molecules"):
            if data.get(key) is not None:
                problems.add(key, "only applies together with 'structure' (cleaning up a complex)")

    if data.get("ligands") is None:
        problems.add("ligands", "missing. Give the SDF file with the ligands to dock")
    else:
        cfg["ligands"] = _path(data["ligands"], "ligands", base, problems, (".sdf",))
    if data.get("output") is None:
        problems.add("output", "missing. Give the folder the results are written to")
    else:
        cfg["output"] = _path(data["output"], "output", base, problems, must_exist=False)

    if data.get("existing_output") is not None:
        cfg["existing_output"] = _choice(data["existing_output"], "existing_output", ("resume", "overwrite"), problems)
    if data.get("protonation") is not None:
        cfg["protonation"] = _choice(data["protonation"], "protonation", ("pdb2pqr", "adfr", "none"), problems)
    if data.get("ph") is not None:
        cfg["ph"] = _number(data["ph"], "ph", problems, float, 0.0, 14.0)
        if cfg.get("protonation") not in (None, "pdb2pqr"):
            problems.add("ph", "only applies with 'protonation: pdb2pqr'")
    if data.get("flagged_ligands") is not None:
        cfg["flagged_ligands"] = _choice(data["flagged_ligands"], "flagged_ligands", ("exclude", "include"), problems)

    cfg["docking"] = _read_docking(_section(data, "docking", DOCKING_KEYS, problems), problems)
    cfg["run"] = _read_run(_section(data, "run", RUN_KEYS, problems), base, problems)

    if problems.lines:
        raise RunFileError(f"Problems in {path}:\n" + "\n".join(f"  - {line}" for line in problems.lines))
    return cfg


def _read_docking(section, problems):
    out = {}
    if "anchor_mode" in section:
        out["anchor_mode"] = _choice(section["anchor_mode"], "docking.anchor_mode", ("spread", "random", "area"), problems)
    if section.get("max_anchors") not in (None, "auto"):
        out["max_anchors"] = _number(section["max_anchors"], "docking.max_anchors", problems, int, 1)
    if "box_size" in section:
        out["box_size"] = _number(section["box_size"], "docking.box_size", problems, float, 1.0)
    if "samples" in section:
        out["samples"] = _number(section["samples"], "docking.samples", problems, int, 1)
    if section.get("seed") is not None:
        out["seed"] = _number(section["seed"], "docking.seed", problems, int, 0)
    if section.get("core_smarts") is not None:
        out["core_smarts"] = _text(section["core_smarts"], "docking.core_smarts", problems)
    if section.get("anchor_atoms") is not None:
        atoms = section["anchor_atoms"]
        if isinstance(atoms, list) and all(isinstance(a, int) and not isinstance(a, bool) and a >= 0 for a in atoms):
            out["anchor_atoms"] = atoms
        else:
            problems.add("docking.anchor_atoms", "must be a list of atom positions (0-based), for example [0, 3, 5]")
        if "core_smarts" not in out:
            problems.add("docking.anchor_atoms", "only applies together with 'core_smarts'")
    if section.get("clean_intermediates") not in (None, "auto"):
        out["clean_intermediates"] = _flag(section["clean_intermediates"], "docking.clean_intermediates", problems)
    if "verbose" in section:
        out["verbose"] = _flag(section["verbose"], "docking.verbose", problems)
    return out


def _read_run(section, base, problems):
    out = {}
    if section.get("mode") is None:
        problems.add("run.mode", "missing. Say what to do with the job: 'submit' (send it to the cluster queue), "
                                 "'run' (run it now in this terminal, on a machine with a GPU) or 'save' (only "
                                 "write the job script)")
    else:
        out["mode"] = _choice(section["mode"], "run.mode", ("submit", "run", "save"), problems)
    if section.get("job_name") is not None:
        name = _text(section["job_name"], "run.job_name", problems)
        if name and not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
            problems.add("run.job_name", "use only letters, digits, '_', '.' and '-'")
        elif name:
            out["job_name"] = name
    if section.get("time") is not None:
        time_limit = section["time"]
        if not isinstance(time_limit, str):
            # YAML reads an unquoted 2:00:00 as the number 7200, which SLURM would take as 7200 minutes
            problems.add("run.time", f"put the time in quotes, for example \"2:00:00\" (without quotes YAML turns it "
                                     f"into the number {time_limit})")
        elif re.fullmatch(r"(\d+-)?\d+(:\d{2}){0,2}", time_limit.strip()):
            out["time"] = time_limit.strip()
        else:
            problems.add("run.time", "use hours:minutes:seconds in quotes, for example \"2:00:00\" "
                                     "(or days-hours:minutes:seconds)")
    if section.get("partition") is not None:
        out["partition"] = _text(section["partition"], "run.partition", problems)
    if section.get("reservation") is not None:
        reservation = _text(section["reservation"], "run.reservation", problems)
        out["reservation"] = "" if reservation and reservation.lower() == "none" else reservation
    if section.get("gpus") is not None:
        out["gpus"] = _number(section["gpus"], "run.gpus", problems, int, 1)
    if section.get("modules") is not None:
        modules = section["modules"]
        if isinstance(modules, list):
            modules = ",".join(str(m) for m in modules)
        out["modules"] = _text(modules, "run.modules", problems)
    if section.get("pixi_dir") is not None:
        out["pixi_dir"] = _path(section["pixi_dir"], "run.pixi_dir", base, problems, must_exist=False)
    if section.get("autodock_gpu") is not None:
        out["autodock_gpu"] = _choice(section["autodock_gpu"], "run.autodock_gpu", ("auto", "cluster"), problems)
    if "watch" in section:
        out["watch"] = _flag(section["watch"], "run.watch", problems)
    return out


# --------------------------------------------------------------------------
# Writing
# --------------------------------------------------------------------------


def _scalar(value):
    text = yaml.safe_dump(value, default_flow_style=True, width=10 ** 6).strip()
    return text[:-4].rstrip() if text.endswith("\n...") else text


def _line(key, value, comment="", indent=0):
    text = f"{' ' * indent}{key}: {_scalar(value)}"
    return f"{text.ljust(44)} # {comment}" if comment else text


def render(values, run_command="pixi run dock"):
    """Write the answers of a wizard run as a run file. `values` has the keys of the file (see load()), with the
    sections 'docking' and 'run' as dicts. Keys with the value None are left out."""
    lines = [f"# Anchored docking run file, written by the wizard. Run it again with:  {run_command} <this file>",
             "# Edit the values to change the run; paths without a leading / are relative to this file.", ""]
    inputs = [("structure", "complex to clean up (PDB or mmCIF)"),
              ("template_molecule", "the molecule to anchor to: residue code, chain and number"),
              ("chains", "protein chains kept"),
              ("keep_molecules", "ions or cofactors kept in the receptor ([] for none)"),
              ("receptor", "protein PDB file"),
              ("template", "template ligand SDF file")]
    for key, comment in inputs:
        if values.get(key) is not None:
            lines.append(_line(key, values[key], comment))
    lines.append(_line("ligands", values["ligands"], "ligands to dock (SDF)"))
    lines.append(_line("output", values["output"], "results folder"))
    for key, comment in (("existing_output", "output folder is not empty: resume | overwrite"),
                         ("protonation", "receptor had no hydrogens: pdb2pqr | adfr | none"),
                         ("ph", "pH for pdb2pqr"),
                         ("flagged_ligands", "ligands that were flagged: exclude | include")):
        if values.get(key) is not None:
            lines.append(_line(key, values[key], comment))
    for name, comments in (("docking", {
            "anchor_mode": "spread | random | area", "max_anchors": "most anchor atoms per ligand",
            "box_size": "search box in Angstrom", "samples": "random attempts ('area' mode)", "seed": "random seed",
            "core_smarts": "manual core pattern", "anchor_atoms": "anchor positions in that pattern (0-based)",
            "clean_intermediates": "delete grid maps etc. when a ligand is done", "verbose": "extra debug output"}),
            ("run", {
            "mode": "submit | run | save", "job_name": "", "time": "time limit, hours:minutes:seconds",
            "partition": "", "reservation": "", "gpus": "", "modules": "loaded in the job, GCC first",
            "pixi_dir": "folder with the pixi program", "autodock_gpu": "auto | cluster", "watch": "wait for the job"})):
        section = {k: v for k, v in (values.get(name) or {}).items() if v is not None}
        if section:
            lines += ["", f"{name}:"]
            lines += [_line(k, v, comments.get(k, ""), 2) for k, v in section.items()]
    return "\n".join(lines) + "\n"


EXAMPLE = """\
# Anchored docking run file. Start it with:  dock run.yaml      (check it first with:  dock run.yaml --check)
# Paths without a leading / are relative to this file. Everything the wizard would ask is either set here or has
# the default noted below. Where the data forces a decision that is missing, the run stops and says what to add.

# --- 1. The protein: EITHER a complex to clean up ... ---------------------------------------------------------
structure: 4djp.pdb              # PDB or mmCIF file of a protein-ligand complex
template_molecule: B7Z A401      # the molecule to anchor to: residue code, chain and number (see the error message
                                 # for the candidates if you do not know them)
chains: [A]                      # protein chains to keep (needed when the file has more than one)
keep_molecules: []               # ions or cofactors near the template to keep, e.g. [ZN A601]; [] keeps none

# --- ... OR a receptor and template that are ready to use -----------------------------------------------------
# receptor: protein.pdb          # protein PDB file, without the ligand
# template: template.sdf         # template ligand SDF file, positioned in the receptor

# --- 2. What to dock and where the results go -----------------------------------------------------------------
ligands: ligands.sdf             # SDF file with the ligands that share a core with the template
output: results/run1             # results folder
# existing_output: resume        # needed only if the folder is not empty: resume (skip finished ligands) | overwrite

# --- 3. Decisions the wizard asks about, needed only when they arise ------------------------------------------
# protonation: pdb2pqr           # needed if the receptor has no hydrogens: pdb2pqr (predicts charges at a pH) | adfr (quick guess) | none
# ph: 7.4                        # pH for pdb2pqr (default 7.4)
# flagged_ligands: exclude       # needed if some ligands are flagged: exclude | include

# --- 4. Docking settings (all optional; these are the wizard's recommendations) -------------------------------
docking:
  anchor_mode: spread            # spread | random | area
  # max_anchors: 32              # default: as many as AutoDock-GPU supports
  box_size: 20                   # search box in Angstrom
  # samples: 100                 # random attempts for 'area' mode
  # seed: 0                      # random seed, for 'random' and 'area'
  # core_smarts: "c1ccccc1"      # manual core pattern, instead of finding the core automatically
  # anchor_atoms: [0, 1, 2]      # which atoms of that pattern are anchors (0-based)
  # clean_intermediates: true    # default: true for 10 or more ligands
  # verbose: false

# --- 5. How to run it -----------------------------------------------------------------------------------------
run:
  mode: submit                   # submit (cluster queue) | run (here, needs a GPU) | save (only write the script)
  # job_name: dock_ligands       # default: from the ligand file name
  # time: "2:00:00"              # default: 10 minutes per ligand, at least 30 minutes
  # partition: ga100
  # reservation: none
  # gpus: 1
  # modules: GCC/13.2.0,AutoDock-GPU/1.5.3-CUDA
  # pixi_dir: /home/me/.local/bin    # default: where pixi is found
  # autodock_gpu: auto           # auto | cluster (the cluster's module, even if a ligand then gets fewer anchors)
  # watch: false                 # true: wait for the job to finish and show the results
"""
