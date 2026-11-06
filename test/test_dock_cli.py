"""Tests for the interactive CLI (dock_cli.py).

Run with:  pixi run -e dev test

None of these need a GPU, SLURM or AutoDock-GPU. SLURM calls are mocked.
"""

import builtins
import json
import os
import shutil
import stat
import subprocess
import sys

import pytest
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

import anchored_docking as ad
import dock_cli as d

RDLogger.DisableLog("rdApp.*")

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, "test")
RECEPTOR = os.path.join(DATA, "Bace1_protein.pdb")
TEMPLATE = os.path.join(DATA, "crystal_ligand.sdf")
LIGAND = os.path.join(DATA, "CAT-13f.sdf")


@pytest.fixture(autouse=True)
def clear_warnings():
    d.WARNINGS.clear()
    yield
    d.WARNINGS.clear()


@pytest.fixture(scope="module")
def template():
    return d.load_template(TEMPLATE)


@pytest.fixture(scope="module")
def receptor():
    return d.analyse_receptor(RECEPTOR)


@pytest.fixture(scope="module")
def reference_ligand():
    return Chem.SDMolSupplier(LIGAND)[0]


def write_sdf(path, mols):
    writer = Chem.SDWriter(str(path))
    writer.SetKekulize(False)
    for mol in mols:
        writer.write(mol)
    writer.close()
    return str(path)


def named(smiles, name):
    mol = Chem.MolFromSmiles(smiles)
    mol.SetProp("_Name", name)
    return mol


def warnings_text():
    return "\n".join(d.WARNINGS)


def atom_line(record, serial, name, resname, x, y, z, element, altloc=" ", chain="A", resseq=1):
    return (f"{record:<6}{serial:5d} {name:<4}{altloc}{resname:>3} {chain}{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}\n")


def edited_receptor(tmp_path, keep=None, extra=()):
    """Copy the example receptor, optionally filtering its lines and appending more."""
    out = tmp_path / "receptor.pdb"
    with open(RECEPTOR) as src, open(out, "w") as dst:
        for line in src:
            if keep is None or keep(line):
                dst.write(line)
        dst.writelines(extra)
    return str(out)


# --------------------------------------------------------------------------
# Receptor checks
# --------------------------------------------------------------------------


def test_example_receptor_is_clean(receptor, capsys):
    assert d.report_receptor(receptor) is True
    assert d.WARNINGS == []
    assert receptor["hydrogens"] > 0
    assert "heavy atoms" in capsys.readouterr().out


def test_receptor_without_hydrogens_warns(tmp_path):
    path = edited_receptor(tmp_path, keep=lambda l: not (l.startswith("ATOM") and l[76:78].strip() == "H"))
    assert d.report_receptor(d.analyse_receptor(path)) is True
    assert "NO hydrogen" in warnings_text()


def test_receptor_with_ligand_and_water_warns_only_about_the_ligand(tmp_path):
    extra = [atom_line("HETATM", 9001, "C1", "LIG", 1.0, 2.0, 3.0, "C"),
             atom_line("HETATM", 9002, "O", "HOH", 4.0, 5.0, 6.0, "O")]
    info = d.analyse_receptor(edited_receptor(tmp_path, extra=extra))
    assert info["waters"] == 1 and info["hetnames"] == {"LIG"}
    d.report_receptor(info)
    assert "LIG" in warnings_text() and "HOH" not in warnings_text()


def test_receptor_with_several_models_warns(tmp_path):
    path = edited_receptor(tmp_path, extra=["MODEL        1\n", "MODEL        2\n"])
    d.report_receptor(d.analyse_receptor(path))
    assert "2 models" in warnings_text()


def test_receptor_with_alternate_locations_warns(tmp_path):
    extra = [atom_line("ATOM", 9001, "CA", "GLY", 1.0, 2.0, 3.0, "C", altloc="B")]
    d.report_receptor(d.analyse_receptor(edited_receptor(tmp_path, extra=extra)))
    assert "alternate locations" in warnings_text()


def test_file_without_atoms_is_rejected(tmp_path, capsys):
    path = tmp_path / "notes.pdb"
    path.write_text("this is not a structure\n")
    assert d.report_receptor(d.analyse_receptor(str(path))) is False
    assert "No protein atoms" in capsys.readouterr().out


# --------------------------------------------------------------------------
# Template checks
# --------------------------------------------------------------------------


def moved(mol, dx):
    mol = Chem.Mol(mol)
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x + dx, p.y, p.z))
    return mol


def test_example_template_sits_in_receptor(template, receptor, capsys):
    extent = d.report_template(template, receptor["coords"])
    assert extent == pytest.approx(11.8, abs=0.1)
    assert d.WARNINGS == []
    assert "sits in the receptor" in capsys.readouterr().out


def test_template_far_from_receptor_warns_about_coordinate_systems(template, receptor):
    assert d.report_template(moved(template, 60.0), receptor["coords"]) is not None
    assert "different coordinate systems" in warnings_text()


def test_template_on_top_of_receptor_atoms_warns(template):
    coords = [tuple(p) for p in template.GetConformer().GetPositions()]
    d.report_template(template, coords)
    assert "on top of receptor atoms" in warnings_text()


def test_template_without_3d_coordinates_is_rejected(capsys):
    mol = Chem.MolFromSmiles("CCO")
    AllChem.Compute2DCoords(mol)
    assert d.report_template(mol, []) is None
    assert "3D" in capsys.readouterr().out


def test_flat_coordinates_flagged_as_3d_are_rejected(capsys):
    mol = Chem.MolFromSmiles("CCOC")
    AllChem.Compute2DCoords(mol)
    mol.GetConformer().Set3D(True)
    assert d.report_template(mol, []) is None
    assert "2D" in capsys.readouterr().out


def test_unreadable_template_is_rejected(capsys):
    assert d.report_template(None, []) is None
    assert "could not read" in capsys.readouterr().out


def test_template_file_with_several_molecules_warns(template, receptor, tmp_path):
    path = write_sdf(tmp_path / "two.sdf", [template, template])
    assert d.count_molecules(path) == 2
    d.report_template(template, receptor["coords"], d.count_molecules(path))
    assert "only the first is used" in warnings_text()


# --------------------------------------------------------------------------
# Ligand checks
# --------------------------------------------------------------------------


def test_example_ligand_is_usable_without_notes(template):
    (result,) = d.check_ligands(template, LIGAND)
    assert result.usable and result.notes == []
    assert result.anchors >= 3 and result.core_atoms >= d.MIN_CORE_ATOMS
    assert result.name == "CAT-13f"


def test_ligand_from_another_series_is_rejected(template, tmp_path):
    (result,) = d.check_ligands(template, write_sdf(tmp_path / "l.sdf", [named("c1ccccc1CCN", "amine")]))
    assert not result.usable and "no shared core" in result.reason


def test_tiny_shared_core_is_flagged_but_usable(template, tmp_path):
    (result,) = d.check_ligands(template, write_sdf(tmp_path / "l.sdf", [named("CCO", "ethanol")]))
    assert result.usable and any("shares only" in n for n in result.notes)
    assert d.flagged_ligands([result]) == [result]


def test_salt_is_flagged(template, reference_ligand, tmp_path):
    salt = Chem.CombineMols(reference_ligand, Chem.MolFromSmiles("[Cl-]"))
    salt.SetProp("_Name", "salt")
    (result,) = d.check_ligands(template, write_sdf(tmp_path / "l.sdf", [salt]))
    assert result.usable and any("separate fragments" in n for n in result.notes)


def test_unreadable_molecule_is_reported_not_fatal(template, tmp_path):
    bad = Chem.MolFromSmiles("C(C)(C)(C)(C)C", sanitize=False)
    bad.SetProp("_Name", "bad")
    path = write_sdf(tmp_path / "l.sdf", [bad, Chem.SDMolSupplier(LIGAND)[0]])
    bad_result, good_result = d.check_ligands(template, path)
    assert not bad_result.usable and "could not read" in bad_result.reason
    assert good_result.usable


def test_atoms_without_an_autodock_type_are_found():
    types = d.load_atom_types()
    assert d.untyped_symbols(Chem.MolFromSmiles("c1ccccc1[Se]C"), types) == ["Se"]
    assert d.untyped_symbols(Chem.MolFromSmiles("CC(=O)Nc1ccc(Cl)cc1"), types) == []
    notes = d.ligand_notes(Chem.MolFromSmiles("c1ccccc1[Se]C"), 12, types)
    assert any("no atom type" in n and "Se" in n for n in notes)


def test_ligand_report_lists_problems(template, tmp_path, capsys):
    mols = [Chem.SDMolSupplier(LIGAND)[0], named("c1ccccc1CCN", "amine"), named("CCO", "ethanol")]
    results = d.check_ligands(template, write_sdf(tmp_path / "l.sdf", mols))
    d.report_ligands(results)
    out = capsys.readouterr().out
    assert "3 ligand(s) in the file; 2 share a core" in out
    assert "amine" in out and "ethanol" in out
    assert [r.name for r in d.flagged_ligands(results)] == ["amine", "ethanol"]


# --------------------------------------------------------------------------
# Limits and prompts
# --------------------------------------------------------------------------


def scripted_input(monkeypatch, answers):
    answers = iter(answers)
    monkeypatch.setattr(builtins, "input", lambda prompt="": next(answers))


def test_box_limit_matches_the_grid_point_limit():
    assert d.grid_points(d.MAX_BOX_SIZE) <= d.MAX_GRID_POINTS
    assert d.grid_points(d.MAX_BOX_SIZE + 1) > d.MAX_GRID_POINTS


def test_anchor_limit_follows_the_atom_type_limit():
    with open(os.path.join(REPO, "scripts", "write-gpf.py")) as fh:
        line = next(l for l in fh if l.startswith("ligand_types "))
    assert d.standard_ligand_types() == len(line.split()) - 1
    assert d.max_anchor_count() == d.MAX_ATOM_TYPES - d.standard_ligand_types()
    assert d.max_anchor_count() >= 3


def test_disk_estimate_grows_with_the_number_of_anchors():
    assert d.estimate_disk_mb(10, 20.0, False, 17) > d.estimate_disk_mb(10, 20.0, False, 3)
    assert d.estimate_disk_mb(10, 20.0, False) == d.estimate_disk_mb(10, 20.0, False, 3)
    # measured on the Bace1 example: 24 MB with 3 anchors, 41 MB with 17 (the estimate is a bit generous)
    assert 30 <= d.disk_per_ligand_mb(20.0, 3) <= 40 and 41 <= d.disk_per_ligand_mb(20.0, 17) <= 60


def test_anchors_in_use_follow_the_mode_and_the_smallest_core():
    class Ligand:
        def __init__(self, anchors):
            self.anchors = anchors

    spread = d.Settings(anchor_mode="spread", anchor_num=17, kept=[Ligand(12), Ligand(25)])
    assert d.anchors_in_use(spread) == 17                      # the larger core is limited to 17
    assert d.anchors_in_use(d.Settings(anchor_mode="spread", anchor_num=17, kept=[Ligand(6), Ligand(9)])) == 9
    assert d.anchors_in_use(d.Settings(anchor_mode="area", anchor_num=3, kept=[])) == 3
    assert d.anchors_in_use(d.Settings(anchor_mode="random", anchor_num=5, kept=[])) == 5


def test_disk_estimate_scales_with_box_volume_and_cleaning():
    assert d.estimate_disk_mb(10, 40.0, False) == pytest.approx(8 * d.estimate_disk_mb(10, 20.0, False))
    assert d.estimate_disk_mb(100, 20.0, True) < d.estimate_disk_mb(100, 20.0, False)


def test_negative_seed_is_rejected_then_valid_one_accepted(monkeypatch, capsys):
    scripted_input(monkeypatch, ["-1", "abc", "5"])
    assert d.ask_number("Random seed", None, int, 0, required=False) == 5
    assert capsys.readouterr().out.count("[x]") == 2


def test_numbers_above_the_maximum_are_rejected(monkeypatch):
    scripted_input(monkeypatch, [str(d.max_anchor_count() + 1), "4"])
    assert d.ask_number("Anchors", 3, int, 1, maximum=d.max_anchor_count()) == 4


def test_blank_seed_means_none(monkeypatch):
    scripted_input(monkeypatch, [""])
    assert d.ask_number("Random seed", None, int, 0, required=False) is None


@pytest.mark.parametrize("text", ["2:00:00", "90:00", "30", "1-12:00:00"])
def test_valid_time_limits(text):
    assert d.time_validator(text) == text


@pytest.mark.parametrize("text", ["2x", "", "1:2:3:4", "two hours", "-5"])
def test_invalid_time_limits(text):
    with pytest.raises(ValueError):
        d.time_validator(text)


def test_minutes_are_formatted_for_slurm():
    assert d.fmt_minutes(30) == "0:30:00"
    assert d.fmt_minutes(200) == "3:20:00"
    assert d.fmt_minutes(1600) == "1-02:40:00"


def test_file_validator(tmp_path):
    ok_file = tmp_path / "a.sdf"
    ok_file.write_text("x")
    empty = tmp_path / "b.sdf"
    empty.write_text("")
    check = d.file_validator((".sdf",))
    assert check(str(ok_file)) == str(ok_file)
    for bad in (str(tmp_path / "missing.sdf"), str(empty), str(tmp_path)):
        with pytest.raises(ValueError):
            check(bad)
    (tmp_path / "c.pdb").write_text("x")
    with pytest.raises(ValueError):
        check(str(tmp_path / "c.pdb"))


# --------------------------------------------------------------------------
# Output folder and environment
# --------------------------------------------------------------------------


def test_output_folder_must_not_contain_inputs_or_repo(tmp_path):
    inputs = tmp_path / "project"
    inputs.mkdir()
    receptor = inputs / "receptor.pdb"
    receptor.write_text("x")
    protected = [str(receptor), REPO, os.path.expanduser("~"), "/"]
    for unsafe in (str(inputs), str(tmp_path), REPO, os.path.expanduser("~"), "/"):
        assert d.validate_outdir(unsafe, protected), unsafe
    assert d.validate_outdir(str(inputs / "output"), protected) is None


def test_output_folder_rejects_spaces_and_files(tmp_path):
    assert "spaces" in d.validate_outdir(str(tmp_path / "my output"), [])
    afile = tmp_path / "file"
    afile.write_text("x")
    assert "file" in d.validate_outdir(str(afile), [])


@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores permissions")
def test_output_folder_must_be_writable(tmp_path):
    locked = tmp_path / "locked"
    locked.mkdir()
    locked.chmod(stat.S_IRUSR | stat.S_IXUSR)
    try:
        if os.access(locked, os.W_OK):
            pytest.skip("file system ignores permission bits")
        assert "permission" in d.validate_outdir(str(locked / "output"), [])
    finally:
        locked.chmod(stat.S_IRWXU)


def test_adfr_detection(tmp_path):
    assert d.adfr_ready(str(tmp_path)) is False
    exe = tmp_path / "bin" / "prepare_receptor"
    exe.parent.mkdir()
    exe.write_text("#!/bin/sh\n")
    assert d.adfr_ready(str(tmp_path)) is False  # present but not executable
    exe.chmod(0o755)
    assert d.adfr_ready(str(tmp_path)) is True


def test_missing_adfr_stops_the_cli_with_instructions(monkeypatch, capsys):
    monkeypatch.setattr(d, "adfr_ready", lambda path=None: False)
    with pytest.raises(SystemExit) as info:
        d.check_environment()
    assert info.value.code == 1
    assert "pixi install --all" in capsys.readouterr().out


# --------------------------------------------------------------------------
# Manual core pattern
# --------------------------------------------------------------------------


def template_scaffold(template):
    from rdkit.Chem.Scaffolds import MurckoScaffold

    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(template))


def test_manual_core_accepts_a_fragment_of_the_template(template):
    pattern = template_scaffold(template)
    anchors = d.manual_core_anchors(pattern, template)
    assert anchors and all(0 <= a < Chem.MolFromSmiles(pattern).GetNumAtoms() for a in anchors)


def test_automatic_core_is_not_necessarily_a_valid_manual_pattern(template, reference_ligand):
    """The automatic search ignores element types, so its core may not match the template.

    anchored_docking.py would fail with a KeyError in that case; the CLI must say so up front.
    """
    core = ad.find_anchor_points(template, reference_ligand, 0, 100, None).pattern_mol
    pattern = Chem.MolToSmiles(core)
    if template.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
        pytest.skip("for this pair the automatic core happens to match the template")
    with pytest.raises(ValueError, match="does not match"):
        d.manual_core_anchors(pattern, template)


@pytest.mark.parametrize("pattern, message", [("not_a_pattern((", "SMARTS"),
                                              ("[C,N]", "SMILES"),
                                              ("[Pt]", "does not match")])
def test_manual_core_rejects_bad_patterns(template, pattern, message):
    with pytest.raises(ValueError, match=message):
        d.manual_core_anchors(pattern, template)


def test_anchor_positions_are_validated():
    assert d.parse_anchor_indices("1 3", {1, 2, 3}) == [1, 3]
    with pytest.raises(ValueError, match="whole numbers"):
        d.parse_anchor_indices("a b", {1})
    with pytest.raises(ValueError, match="Not usable"):
        d.parse_anchor_indices("1 9", {1, 2})


# --------------------------------------------------------------------------
# SLURM helpers (mocked)
# --------------------------------------------------------------------------


def completed(returncode=0, stdout="", stderr=""):
    return subprocess.CompletedProcess([], returncode, stdout, stderr)


def mock_run(monkeypatch, result):
    def run(*args, **kwargs):
        if isinstance(result, Exception):
            raise result
        return result

    monkeypatch.setattr(d.subprocess, "run", run)


def test_job_state_running(monkeypatch):
    mock_run(monkeypatch, completed(stdout="RUNNING|0:10\n"))
    assert d.job_state("1") == ("state", ["RUNNING", "0:10"])


def test_job_state_gone_when_queue_is_empty_or_job_is_unknown(monkeypatch):
    mock_run(monkeypatch, completed(stdout=""))
    assert d.job_state("1") == ("gone", None)
    mock_run(monkeypatch, completed(1, stderr="slurm_load_jobs error: Invalid job id specified"))
    assert d.job_state("1") == ("gone", None)


def test_squeue_failure_is_an_error_not_a_finished_job(monkeypatch):
    mock_run(monkeypatch, completed(1, stderr="slurm_load_jobs error: Socket timed out"))
    kind, message = d.job_state("1")
    assert kind == "error" and "timed out" in message
    mock_run(monkeypatch, FileNotFoundError("squeue"))
    assert d.job_state("1")[0] == "error"


def test_final_state_is_parsed(monkeypatch):
    mock_run(monkeypatch, completed(stdout="COMPLETED|0:0\n"))
    assert d.job_final_state("1") == ("COMPLETED", "0:0")
    mock_run(monkeypatch, completed(1))
    assert d.job_final_state("1") is None


def watch_settings(tmp_path):
    return d.Settings(outdir=str(tmp_path / "out"), jobname="job")


def test_watch_does_not_claim_success_when_slurm_is_unreachable(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(d.time, "sleep", lambda s: None)
    monkeypatch.setattr(d, "job_state", lambda job: ("error", "Socket timed out"))
    assert d.watch("7", watch_settings(tmp_path)) is False
    out = capsys.readouterr().out
    assert "Stopped watching" in out and "finished" not in out.lower()


def test_watch_reports_a_finished_job(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(d.time, "sleep", lambda s: None)
    states = iter([("state", ["PENDING", "0:00"]), ("state", ["RUNNING", "0:05"]), ("gone", None)])
    monkeypatch.setattr(d, "job_state", lambda job: next(states))
    monkeypatch.setattr(d, "job_final_state", lambda job: ("COMPLETED", "0:0"))
    assert d.watch("7", watch_settings(tmp_path)) is True
    out = capsys.readouterr().out
    assert "waiting in the queue" in out and "running" in out and "COMPLETED" in out


def test_watch_keeps_waiting_while_accounting_still_shows_the_job_active(monkeypatch, tmp_path):
    monkeypatch.setattr(d.time, "sleep", lambda s: None)
    monkeypatch.setattr(d, "job_state", lambda job: ("gone", None))
    finals = iter([("RUNNING", "0:0"), ("COMPLETED", "0:0")])
    monkeypatch.setattr(d, "job_final_state", lambda job: next(finals))
    assert d.watch("7", watch_settings(tmp_path)) is True


def test_watch_survives_a_few_transient_errors(monkeypatch, tmp_path):
    monkeypatch.setattr(d.time, "sleep", lambda s: None)
    states = iter([("error", "x"), ("error", "x"), ("gone", None)])
    monkeypatch.setattr(d, "job_state", lambda job: next(states))
    monkeypatch.setattr(d, "job_final_state", lambda job: ("COMPLETED", "0:0"))
    assert d.watch("7", watch_settings(tmp_path)) is True


# --------------------------------------------------------------------------
# Command and job script
# --------------------------------------------------------------------------


def make_settings(tmp_path, **overrides):
    settings = d.Settings(receptor="/data/my receptor.pdb", template="/data/template.sdf",
                          ligands="/data/ligands.sdf", outdir=str(tmp_path / "out"), pixi_dir="/opt/bin",
                          jobname="myjob", time_limit="1:00:00")
    for key, value in overrides.items():
        setattr(settings, key, value)
    return settings


def test_command_arguments(tmp_path):
    args = d.build_args(make_settings(tmp_path))
    assert args[args.index("--input_file") + 1] == "/data/ligands.sdf"
    assert "--force" not in args and "--clean" not in args and "--random_seed" not in args
    args = d.build_args(make_settings(tmp_path, force=True, clean=True, seed=0, anchor_indices=[1, 2],
                                      query_smarts="c1ccccc1"))
    assert {"--force", "--clean"} <= set(args)
    assert args[args.index("--random_seed") + 1] == "0"  # seed 0 must not be dropped
    assert args[args.index("--anchor_indices") + 1:args.index("--anchor_indices") + 3] == ["1", "2"]


def test_job_script_quotes_paths_and_is_valid_bash(tmp_path):
    script = d.build_script(make_settings(tmp_path))
    assert "'/data/my receptor.pdb'" in script
    assert "module load GCC/13.2.0\nmodule load AutoDock-GPU/1.5.3-CUDA" in script
    assert f"cd {REPO}" in script and "#SBATCH --time=1:00:00" in script
    assert script.index("module load GCC") < script.index("module load AutoDock")
    path = tmp_path / "job.sh"
    path.write_text(script)
    assert subprocess.run(["bash", "-n", str(path)], capture_output=True).returncode == 0


# --------------------------------------------------------------------------
# Embedding check, resume and results display
# --------------------------------------------------------------------------


def test_ligand_that_cannot_be_embedded_is_rejected_before_submission(template, tmp_path, monkeypatch):
    def refuse(mol):
        raise ad.LigandError("RDKit could not generate 3D coordinates for this ligand")

    monkeypatch.setattr(ad, "embed_ligand", refuse)
    (result,) = d.check_ligands(template, LIGAND)
    assert not result.usable and "3D coordinates" in result.reason
    assert d.flagged_ligands([result]) == [result]


def old_run(tmp_path, done=1, **changes):
    """A folder holding an earlier run of the example, with `done` ligands finished."""
    out = tmp_path / "output"
    out.mkdir()
    values = dict(template_file=TEMPLATE, input_file=LIGAND, receptor_file=RECEPTOR,
                  anchor_num=d.max_anchor_count(),
                  anchor_mode="spread", size=20.0, samples=100, random_seed=None, query_smarts=None,
                  anchor_indices=None)
    values.update(changes)
    (out / "run_settings.json").write_text(json.dumps(ad.run_signature(**values)))
    for i in range(done):
        (out / str(i)).mkdir()
        (out / str(i) / "ligand_docked.sdf").write_text("x")
    return out


def test_wizard_offers_to_resume_a_matching_earlier_run(monkeypatch, tmp_path, capsys):
    out = old_run(tmp_path, done=1)
    run_wizard(monkeypatch, [RECEPTOR, TEMPLATE, LIGAND, "", "", save_choice(), str(out),
                             "",                       # resume: yes
                             "", "", "", ""])          # job name, time limit, change settings, confirm
    text = capsys.readouterr().out
    assert "same inputs and settings: 1 of 1" in text and "resuming: 1 finished" in text
    script = (tmp_path / "submit_dock_CAT-13f.sh").read_text()
    assert "--resume" in script and "--force" not in script
    assert (out / "0" / "ligand_docked.sdf").read_text() == "x"  # nothing was deleted


def test_wizard_will_not_resume_when_settings_changed(monkeypatch, tmp_path, capsys):
    out = old_run(tmp_path, done=1, anchor_num=5)
    run_wizard(monkeypatch, [RECEPTOR, TEMPLATE, LIGAND, "", "", save_choice(), str(out),
                             "yes",                    # no resume is offered, so this answers 'delete contents?'
                             "", "", "", ""])
    text = capsys.readouterr().out
    assert "cannot be resumed: anchor_num" in text
    script = (tmp_path / "submit_dock_CAT-13f.sh").read_text()
    assert "--force" in script and "--resume" not in script


def test_resume_can_be_declined_in_favour_of_starting_again(monkeypatch, tmp_path):
    out = old_run(tmp_path, done=1)
    run_wizard(monkeypatch, [RECEPTOR, TEMPLATE, LIGAND, "", "", save_choice(), str(out),
                             "no",                     # do not resume
                             "yes",                    # delete the contents instead
                             "", "", "", ""])
    script = (tmp_path / "submit_dock_CAT-13f.sh").read_text()
    assert "--force" in script and "--resume" not in script


def test_command_has_resume_but_never_force_together(tmp_path):
    args = d.build_args(make_settings(tmp_path, force=True, resume=True))
    assert "--resume" in args and "--force" not in args


def results_settings(tmp_path, names=("a", "b", "c", "d")):
    kept = [d.LigandResult(i, n, None, True) for i, n in enumerate(names)]
    out = tmp_path / "out"
    out.mkdir()
    return d.Settings(outdir=str(out), kept=kept, all_count=len(kept), jobname="job"), out


def test_results_table_shows_flagged_and_failed_ligands(tmp_path, capsys):
    s, out = results_settings(tmp_path)
    ad.write_summary(str(out), {
        0: {"index": 0, "name": "a", "status": "docked", "anchor_max_dev": 0.4, "core_rmsd": 0.9, "note": ""},
        1: {"index": 1, "name": "b", "status": "flagged", "anchor_max_dev": 4.2, "core_rmsd": 5.0,
            "note": "an anchor atom is 4.2 A from its target"},
        2: {"index": 2, "name": "c", "status": "failed", "note": "autodock_gpu_128wi failed (exit code 1)"}})
    (out / "0").mkdir()
    (out / "0" / "ligand_docked.sdf").write_text("x")
    (out / "1").mkdir()
    (out / "1" / "ligand_docked.sdf").write_text("x")
    d.show_results(s)
    text = capsys.readouterr().out
    assert "CHECK POSE" in text and "FAILED" in text and "not run" in text
    assert "2 of 4 ligand(s) have a docked pose" in text
    assert "1 pose(s) need a closer look" in text and "4.2 A from its target" in text
    assert "1 ligand(s) failed" in text and "resume" in text


def test_results_fall_back_to_files_when_there_is_no_summary(tmp_path, capsys):
    s, out = results_settings(tmp_path, names=("a", "b"))
    (out / "0").mkdir()
    (out / "0" / "ligand_docked.sdf").write_text("x")
    d.show_results(s)
    text = capsys.readouterr().out
    assert "docked (not checked)" in text and "not run" in text


def test_run_here_reports_a_missing_program_without_crashing(monkeypatch, tmp_path):
    def missing(path=None):
        raise FileNotFoundError("AutoDock-GPU binary not found in PATH")

    monkeypatch.setattr(ad, "find_autodock_gpu", missing)
    assert d.run_here(make_settings(tmp_path)) is None


# --------------------------------------------------------------------------
# Whole wizard, driven by scripted answers (save-only mode: nothing is submitted)
# --------------------------------------------------------------------------


def run_wizard(monkeypatch, answers, from_complex=False):
    """Run the whole wizard on scripted answers. By default the first question is answered with
    'separate receptor and template files'; from_complex=True leaves it to the answers."""
    monkeypatch.setattr(sys, "argv", ["dock_cli.py"])
    scripted_input(monkeypatch, answers if from_complex else ["2", *answers])
    d.main()


def save_choice():
    """Menu number of 'only save the script' (SLURM adds a first option when present)."""
    return "3" if shutil.which("sbatch") else "2"


def test_wizard_saves_a_job_script_for_the_example(monkeypatch, tmp_path, capsys):
    out = tmp_path / "output"
    run_wizard(monkeypatch, [RECEPTOR, TEMPLATE, LIGAND, "", "",         # files, recommended settings, clean
                             save_choice(), str(out),                    # how to run, output folder
                             "", "", "",                                 # job name, time limit, change settings
                             ""])                                        # confirm
    script = (tmp_path / "submit_dock_CAT-13f.sh").read_text()
    assert f"--docking_dir {out}" in script
    assert f"--receptor_file {RECEPTOR}" in script
    assert f"--anchor_num {d.max_anchor_count()}" in script and "--anchor_mode spread" in script
    assert "Job script saved" in capsys.readouterr().out
    assert not out.exists()  # saving a script must not run anything


def test_wizard_recommends_spreading_anchors_over_the_whole_core(monkeypatch, tmp_path, capsys):
    run_wizard(monkeypatch, [RECEPTOR, TEMPLATE, LIGAND, "", "", save_choice(), str(tmp_path / "o"), "", "", "", ""])
    text = capsys.readouterr().out
    assert "spread across the shared core" in text and f"at most {d.max_anchor_count()}" in text
    assert "up to" in text and "per ligand, spread across the shared core" in text      # the review table


def test_wizard_asks_for_a_maximum_in_spread_mode_and_no_seed(monkeypatch, tmp_path, capsys):
    run_wizard(monkeypatch, [RECEPTOR, TEMPLATE, LIGAND,
                             "no",                  # not the recommended settings
                             "",                    # anchor selection: spread (default)
                             "99", "5",             # too many, then fine
                             "",                    # box size
                             "",                    # expert options: no (a seed is not asked for in spread mode)
                             "",                    # clean intermediates
                             save_choice(), str(tmp_path / "o"), "", "", "", ""])
    text = capsys.readouterr().out
    assert "uses all of them" in text and "random seed makes" not in text       # the prompt text itself is not captured
    script = (tmp_path / "submit_dock_CAT-13f.sh").read_text()
    assert "--anchor_num 5" in script and "--anchor_mode spread" in script and "--random_seed" not in script


def test_wizard_leaves_out_flagged_ligands_and_validates_expert_input(monkeypatch, tmp_path, capsys):
    mixed = write_sdf(tmp_path / "mixed.sdf", [Chem.SDMolSupplier(LIGAND)[0],
                                               named("c1ccccc1CCN", "amine"), named("CCO", "ethanol")])
    out = tmp_path / "output"
    run_wizard(monkeypatch, [
        str(tmp_path / "missing.pdb"), RECEPTOR, TEMPLATE, mixed,       # first receptor path is wrong
        "",                                                             # leave out flagged ligands: yes
        "no",                                                           # not the recommended settings
        "area", "abc", "50",                                            # anchor mode, bad then good attempts
        "500", "40",                                                    # box too large, then fine
        "-3", "7",                                                      # negative seed, then fine
        "yes",                                                          # expert options
        "[Pt]", "",                                                     # pattern not in the template, then blank
        "",                                                             # no debug output
        "",                                                             # clean intermediates: default
        save_choice(), str(tmp_path),                                   # output folder that holds the inputs
        str(out),                                                       # then a safe one
        "bad name!", "2x", "1:00:00", "",                               # job name, bad then good time, no changes
        ""])                                                            # confirm
    text = capsys.readouterr().out
    assert "File not found" in text
    assert "no shared core" in text and "shares only" in text          # the two flagged ligands
    assert "at most" in text and "does not match the template" in text
    assert "contains" in text and "Choose another folder" in text
    script = (tmp_path / "submit_bad_name_.sh").read_text()
    assert "--anchor_mode area" in script and "--samples 50" in script
    assert "--size 40.0" in script and "--random_seed 7" in script
    reduced = tmp_path / "bad_name__ligands.sdf"
    assert len([m for m in Chem.SDMolSupplier(str(reduced))]) == 1
    assert f"--input_file {reduced}" in script


def test_wizard_can_be_cancelled_at_the_review(monkeypatch, tmp_path, capsys):
    run_wizard(monkeypatch, [RECEPTOR, TEMPLATE, LIGAND, "", "", save_choice(), str(tmp_path / "output"),
                             "", "", "", "no"])
    assert "Cancelled. Nothing was run." in capsys.readouterr().out
    assert list(tmp_path.iterdir()) == []
