"""Tests for the run file: `dock run.yaml` does what the wizard does, with the answers taken from a YAML file.

Run with:  pixi run -e dev test

Nothing here needs a GPU, SLURM or AutoDock-GPU (jobs are only saved, never submitted).
"""

import builtins
import os
import sys

import pytest
from rdkit import Chem

import dock_cli as d
import runfile
import structure_prep as sp
from test_dock_cli import LIGAND, RECEPTOR, TEMPLATE, clear_warnings, run_wizard, save_choice, write_sdf  # noqa: F401
from test_structure_prep import run_complex_wizard, tail, write_complex


def start(monkeypatch, *args):
    """Run `dock_cli.py <args>`; a run file must never ask a question."""
    def no_questions(prompt=""):
        raise AssertionError(f"a run file must not ask questions, but it asked: {prompt}")

    monkeypatch.setattr(sys, "argv", ["dock_cli.py", *args])
    monkeypatch.setattr(builtins, "input", no_questions)
    d.main()


def stops(monkeypatch, capsys, *args):
    """Run a run file that must stop with an explanation. Returns what was printed."""
    with pytest.raises(SystemExit) as stop:
        start(monkeypatch, *args)
    assert stop.value.code == 2
    return capsys.readouterr().out


def write_run(tmp_path, text, name="run.yaml"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def files_run(tmp_path, extra="", mode="save", ligands=LIGAND, receptor=RECEPTOR, output=None):
    """A run file for ready-made receptor and template files."""
    out = output or tmp_path / "work" / "output"
    return write_run(tmp_path, f"receptor: {receptor}\ntemplate: {TEMPLATE}\nligands: {ligands}\noutput: {out}\n"
                               f"{extra}run:\n  mode: {mode}\n")


def without_hydrogens(tmp_path):
    path = tmp_path / "no_h.pdb"
    path.write_text("".join(l for l in open(RECEPTOR) if not (l.startswith("ATOM") and l[76:78].strip() == "H")))
    return str(path)


# --------------------------------------------------------------------------
# Reading the file
# --------------------------------------------------------------------------


def test_the_printed_example_is_a_valid_run_file(tmp_path):
    (tmp_path / "4djp.pdb").write_text("x")
    (tmp_path / "ligands.sdf").write_text("x")
    cfg = runfile.load(write_run(tmp_path, runfile.EXAMPLE))
    assert cfg["route"] == "structure" and cfg["template_molecule"] == "B7Z A401"
    assert cfg["docking"] == {"anchor_mode": "spread", "box_size": 20.0}
    assert cfg["run"] == {"mode": "submit"}


def test_every_problem_is_listed_not_only_the_first(tmp_path):
    path = write_run(tmp_path, f"receptor: {RECEPTOR}\ntemplate: {TEMPLATE}\nligands: missing.sdf\noutput: o\n"
                               "protonation: adfr\nph: 6\nflagged_ligands: maybe\ndocking:\n  anchor_mod: spread\n"
                               "  box_size: -1\nrun:\n  time: soon\n")
    with pytest.raises(runfile.RunFileError) as problem:
        runfile.load(path)
    text = str(problem.value)
    for expected in ("ligands: file not found", "ph: only applies with 'protonation: pdb2pqr'",
                     "flagged_ligands: must be one of: exclude, include",
                     "docking.anchor_mod: unknown key. Did you mean 'anchor_mode'?",
                     "docking.box_size: must be at least 1.0", "run.mode: missing", "run.time:"):
        assert expected in text


def test_a_time_without_quotes_is_refused_because_yaml_turns_it_into_a_number(tmp_path):
    base = f"receptor: {RECEPTOR}\ntemplate: {TEMPLATE}\nligands: {LIGAND}\noutput: o\nrun:\n  mode: save\n"
    with pytest.raises(runfile.RunFileError, match="put the time in quotes.*7200"):
        runfile.load(write_run(tmp_path, base + "  time: 2:00:00\n"))
    assert runfile.load(write_run(tmp_path, base + "  time: \"2:00:00\"\n", "q.yaml"))["run"]["time"] == "2:00:00"


def test_the_time_the_wizard_writes_survives_a_round_trip(tmp_path):
    values = {"receptor": RECEPTOR, "template": TEMPLATE, "ligands": LIGAND, "output": "o",
              "run": {"mode": "save", "time": "2:00:00"}}
    text = runfile.render(values)
    assert "time: '2:00:00'" in text
    assert runfile.load(write_run(tmp_path, text))["run"]["time"] == "2:00:00"


def test_the_protein_must_come_from_one_route(tmp_path):
    both = write_run(tmp_path, f"structure: {RECEPTOR}\nreceptor: {RECEPTOR}\nligands: {LIGAND}\noutput: o\n"
                              "run:\n  mode: save\n", "both.yaml")
    neither = write_run(tmp_path, f"ligands: {LIGAND}\noutput: o\nrun:\n  mode: save\n", "neither.yaml")
    with pytest.raises(runfile.RunFileError, match="not both"):
        runfile.load(both)
    with pytest.raises(runfile.RunFileError, match="say what to dock into"):
        runfile.load(neither)


def test_a_file_that_is_not_yaml_or_not_a_mapping_is_explained(tmp_path):
    with pytest.raises(runfile.RunFileError, match="not valid YAML"):
        runfile.load(write_run(tmp_path, "a: [unclosed\n"))
    with pytest.raises(runfile.RunFileError, match="key: value"):
        runfile.load(write_run(tmp_path, "- just\n- a list\n", "list.yaml"))
    with pytest.raises(runfile.RunFileError, match="Cannot read"):
        runfile.load(str(tmp_path / "absent.yaml"))


def test_paths_are_relative_to_the_run_file_not_to_the_working_directory(tmp_path, monkeypatch):
    folder = tmp_path / "project"
    folder.mkdir()
    (folder / "protein.pdb").write_text("x")
    (folder / "template.sdf").write_text("x")
    (folder / "ligands.sdf").write_text("x")
    path = write_run(folder, "receptor: protein.pdb\ntemplate: template.sdf\nligands: ligands.sdf\noutput: results\n"
                             "run:\n  mode: save\n")
    monkeypatch.chdir(tmp_path)
    cfg = runfile.load(path)
    assert cfg["receptor"] == str(folder / "protein.pdb") and cfg["output"] == str(folder / "results")


def test_what_the_wizard_writes_is_read_back(tmp_path):
    values = {"receptor": RECEPTOR, "template": TEMPLATE, "ligands": LIGAND, "output": str(tmp_path / "o"),
              "existing_output": "overwrite", "protonation": "pdb2pqr", "ph": 7.0, "flagged_ligands": "exclude",
              "docking": {"anchor_mode": "random", "max_anchors": 4, "box_size": 25.0, "seed": 0,
                          "clean_intermediates": False},
              "run": {"mode": "submit", "job_name": "j", "time": "2:00:00", "partition": "ga100",
                      "reservation": "none", "gpus": 1, "modules": "GCC/13.2.0,AutoDock-GPU/1.5.3-CUDA"}}
    cfg = runfile.load(write_run(tmp_path, runfile.render(values, "ddsuite dock")))
    assert cfg["protonation"] == "pdb2pqr" and cfg["ph"] == 7.0 and cfg["flagged_ligands"] == "exclude"
    assert cfg["docking"] == {"anchor_mode": "random", "max_anchors": 4, "box_size": 25.0, "seed": 0,
                              "clean_intermediates": False}
    assert cfg["run"]["reservation"] == "" and cfg["run"]["modules"].startswith("GCC")
    assert "ddsuite dock <this file>" in runfile.render(values, "ddsuite dock")


# --------------------------------------------------------------------------
# Running it
# --------------------------------------------------------------------------


def test_a_run_file_saves_the_job_without_asking_anything(monkeypatch, tmp_path, capsys):
    path = files_run(tmp_path)
    start(monkeypatch, path)
    work = tmp_path / "work"
    script = (work / "submit_dock_CAT-13f.sh").read_text()
    assert f"--docking_dir {work / 'output'}" in script and f"--receptor_file {RECEPTOR}" in script
    assert f"--anchor_num {d.MAX_ATOM_TYPES}" in script and "--anchor_mode spread" in script
    assert not (work / "output").exists()                      # saving a script runs nothing
    saved = work / "dock_CAT-13f_run.yaml"                       # a copy with every default filled in
    assert saved.is_file() and "anchor_mode: spread" in saved.read_text()
    assert "Run file:" in capsys.readouterr().out


def test_check_goes_through_every_check_and_runs_nothing(monkeypatch, tmp_path, capsys):
    start(monkeypatch, files_run(tmp_path), "--check")
    text = capsys.readouterr().out
    assert "Nothing was run (--check)" in text and "Review" in text
    assert not (tmp_path / "work").exists()


def test_check_needs_a_run_file(monkeypatch):
    with pytest.raises(SystemExit) as stop:
        start(monkeypatch, "--check")
    assert stop.value.code == 2


def test_example_prints_a_run_file_and_needs_no_data(monkeypatch, capsys):
    start(monkeypatch, "--example")
    assert capsys.readouterr().out == runfile.EXAMPLE


def test_the_run_file_the_wizard_writes_repeats_the_wizard_job(monkeypatch, tmp_path, capsys):
    out = tmp_path / "output"
    run_wizard(monkeypatch, [RECEPTOR, TEMPLATE, LIGAND, "", "", save_choice(), str(out), "", "", "", ""])
    assert "Run file saved:" in capsys.readouterr().out
    script_path = tmp_path / "submit_dock_CAT-13f.sh"
    run_path = tmp_path / "dock_CAT-13f_run.yaml"
    wizard_script = script_path.read_text()
    assert f"receptor: {RECEPTOR}" in run_path.read_text() and "mode: save" in run_path.read_text()
    script_path.unlink()
    start(monkeypatch, str(run_path))
    assert script_path.read_text() == wizard_script


def test_a_complex_the_wizard_cleaned_is_repeated_from_its_run_file(monkeypatch, tmp_path, capsys):
    structure = write_complex(tmp_path)
    out = tmp_path / "work" / "output"
    run_complex_wizard(monkeypatch, ["1", structure, "", "yes", LIGAND, "", "", *tail(out)])
    work = tmp_path / "work"
    run_path = work / "dock_CAT-13f_run.yaml"
    text = run_path.read_text()
    assert f"structure: {structure}" in text and "template_molecule: LIG" in text and "ZN" in text
    first = {name: (work / name).read_text() for name in
             ("submit_dock_CAT-13f.sh", "dock_CAT-13f_receptor.pdb", "dock_CAT-13f_template.sdf")}
    for name in first:
        (work / name).unlink()
    capsys.readouterr()
    start(monkeypatch, str(run_path))
    assert {name: (work / name).read_text() for name in first} == first


def test_protonation_has_no_default_when_the_receptor_has_no_hydrogens(monkeypatch, tmp_path, capsys):
    path = files_run(tmp_path, receptor=without_hydrogens(tmp_path))
    text = stops(monkeypatch, capsys, path)
    assert "'protonation' is missing" in text and "pdb2pqr" in text and "adfr" in text and "Nothing was run" in text
    assert not (tmp_path / "work").exists()


def test_adfr_protonation_is_accepted(monkeypatch, tmp_path, capsys):
    path = files_run(tmp_path, "protonation: adfr\n", receptor=without_hydrogens(tmp_path))
    start(monkeypatch, path)
    assert (tmp_path / "work" / "submit_dock_CAT-13f.sh").is_file()
    assert "protonation: adfr" in (tmp_path / "work" / "dock_CAT-13f_run.yaml").read_text()


def test_protonation_is_ignored_when_the_receptor_already_has_hydrogens(monkeypatch, tmp_path, capsys):
    start(monkeypatch, files_run(tmp_path, "protonation: adfr\n"))
    assert "is ignored because the receptor already has hydrogens" in capsys.readouterr().out


def test_an_existing_output_folder_must_be_decided(monkeypatch, tmp_path, capsys):
    out = tmp_path / "work" / "output"
    out.mkdir(parents=True)
    (out / "old.txt").write_text("x")
    text = stops(monkeypatch, capsys, files_run(tmp_path))
    assert "'existing_output' is missing" in text and "already exists and is not empty" in text
    assert "resume" in text and "overwrite" in text
    text = stops(monkeypatch, capsys, files_run(tmp_path, "existing_output: resume\n"))
    assert "'existing_output: resume' is not possible" in text
    start(monkeypatch, files_run(tmp_path, "existing_output: overwrite\n"))
    assert "--force" in (tmp_path / "work" / "submit_dock_CAT-13f.sh").read_text()
    assert "existing_output: overwrite" in (tmp_path / "work" / "dock_CAT-13f_run.yaml").read_text()


def test_the_output_folder_may_not_contain_the_inputs(monkeypatch, tmp_path, capsys):
    text = stops(monkeypatch, capsys, files_run(tmp_path, output=os.path.dirname(RECEPTOR)))
    assert "'output'" in text and "Choose another folder" in text


def test_flagged_ligands_must_be_decided(monkeypatch, tmp_path, capsys):
    ligands = tmp_path / "mixed.sdf"
    write_sdf(ligands, [Chem.SDMolSupplier(LIGAND)[0], Chem.MolFromSmiles("CCO")])   # the second shares one atom: flagged
    text = stops(monkeypatch, capsys, files_run(tmp_path, ligands=ligands))
    assert "'flagged_ligands' is missing" in text and "#2 molecule 2: shares only 1 heavy atom" in text
    assert "exclude" in text and "include" in text
    start(monkeypatch, files_run(tmp_path, "flagged_ligands: exclude\n", ligands=ligands))
    work = tmp_path / "work"
    assert "flagged_ligands: exclude" in (work / "dock_mixed_run.yaml").read_text()
    reduced = work / "dock_mixed_ligands.sdf"                      # only the ligand that can be docked
    assert len(list(Chem.SDMolSupplier(str(reduced)))) == 1
    start(monkeypatch, files_run(tmp_path, "flagged_ligands: include\n", ligands=ligands))
    assert "Flagged ligands are included" in capsys.readouterr().out


def test_the_run_mode_cannot_be_left_out(monkeypatch, tmp_path, capsys):
    path = write_run(tmp_path, f"receptor: {RECEPTOR}\ntemplate: {TEMPLATE}\nligands: {LIGAND}\noutput: o\n")
    text = stops(monkeypatch, capsys, path)
    assert "run.mode: missing" in text and "submit" in text and "save" in text and "Nothing was run" in text


def test_a_problem_in_the_file_is_reported_with_the_key(monkeypatch, tmp_path, capsys):
    text = stops(monkeypatch, capsys, files_run(tmp_path, "docking:\n  max_anchors: 99\n"))
    assert "docking.max_anchors is 99" in text and "32" in text


def test_a_complex_run_file_asks_for_the_decisions_the_wizard_asks_for(monkeypatch, tmp_path, capsys):
    structure = write_complex(tmp_path, second_chain=True)
    base = f"structure: {structure}\nligands: {LIGAND}\noutput: {tmp_path / 'work' / 'output'}\n"
    tail_ = "run:\n  mode: save\n"
    text = stops(monkeypatch, capsys, write_run(tmp_path, base + tail_))
    assert "'template_molecule' is missing" in text and "LIG" in text and "the wizard suggests this one" in text
    text = stops(monkeypatch, capsys, write_run(tmp_path, base + "template_molecule: NOPE A1\n" + tail_))
    assert "is not a molecule that can be the template" in text
    inv = sp.inventory(*sp.load_structure(structure), structure)
    label = sp.suggested_template(inv).label
    named = base + f"template_molecule: {label}\n"
    text = stops(monkeypatch, capsys, write_run(tmp_path, named + tail_))
    assert "'chains' is missing" in text and "A" in text and "B" in text
    text = stops(monkeypatch, capsys, write_run(tmp_path, named + "chains: [A]\n" + tail_))
    assert "'keep_molecules' is missing" in text and "ZN" in text
    text = stops(monkeypatch, capsys, write_run(tmp_path, named + "chains: [A]\nkeep_molecules: [XX A9]\n" + tail_))
    assert "'keep_molecules' names XX A9" in text
    start(monkeypatch, write_run(tmp_path, named + "chains: [A]\nkeep_molecules: []\n" + tail_))
    assert (tmp_path / "work" / "submit_dock_CAT-13f.sh").is_file()
