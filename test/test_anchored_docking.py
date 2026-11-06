"""Tests for the batch behaviour of anchored_docking.py: failures, resume and the summary file.

The external programs (AutoGrid, AutoDock-GPU, Meeko, ADFR) are not run: the per-ligand
workflow is replaced by a stub, so these tests need no GPU.
"""

import json
import os
import shutil
import subprocess

import pytest
from rdkit import Chem, RDLogger

import anchored_docking as ad

RDLogger.DisableLog("rdApp.*")

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RECEPTOR = os.path.join(REPO, "test", "Bace1_protein.pdb")
TEMPLATE = os.path.join(REPO, "test", "crystal_ligand.sdf")
LIGAND = os.path.join(REPO, "test", "CAT-13f.sdf")

# Fails 3D embedding with both of RDKit's strategies.
UNEMBEDDABLE = "C1CC[C@]23CC[C@@]1(CC2)CC3"


# --------------------------------------------------------------------------
# Embedding
# --------------------------------------------------------------------------


def test_normal_ligand_is_embedded_with_hydrogens():
    mol = ad.embed_ligand(Chem.SDMolSupplier(LIGAND)[0])
    assert mol.GetNumConformers() == 1
    assert mol.GetNumAtoms() > mol.GetNumHeavyAtoms()


def test_failed_embedding_is_an_error_not_a_flat_ligand():
    with pytest.raises(ad.LigandError, match="3D coordinates"):
        ad.embed_ligand(Chem.MolFromSmiles(UNEMBEDDABLE))


def test_embedding_gives_up_at_the_time_limit():
    """A hopeless ligand must not hold up a job: the time limit is honoured on each attempt."""
    import time

    big_strained = Chem.MolFromSmiles("c1ccccc1CCNC(=O)C1CC[C@]23CC[C@@]1(CC2)CC3")
    started = time.time()
    with pytest.raises(ad.LigandError):
        ad.embed_ligand(big_strained, timeout=1)
    assert time.time() - started < 10  # two attempts of at most about 1 s each, plus overhead


def test_embedding_time_limit_does_not_affect_ordinary_ligands():
    peptide = Chem.MolFromSmiles("CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](CCCN=C(N)N)NC(=O)CN)"
                                 "C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O")
    assert ad.embed_ligand(peptide, timeout=5).GetNumConformers() == 1


def test_failed_embedding_writes_no_ligand_file(tmp_path, monkeypatch):
    called = []
    monkeypatch.setattr(ad.subprocess, "run", lambda *a, **k: called.append(a))
    with pytest.raises(ad.LigandError):
        ad.prepare_ligand(str(tmp_path), Chem.MolFromSmiles(UNEMBEDDABLE))
    assert not (tmp_path / "ligand.sdf").exists() and called == []  # Meeko is never reached


# --------------------------------------------------------------------------
# Run signature and resume rules
# --------------------------------------------------------------------------


def signature(**changes):
    values = dict(template_file=TEMPLATE, input_file=LIGAND, receptor_file=RECEPTOR, anchor_num=3,
                  anchor_mode="random", size=20.0, samples=100, random_seed=None, query_smarts=None,
                  anchor_indices=None)
    values.update(changes)
    return ad.run_signature(**values)


def test_signature_only_changes_with_relevant_inputs(tmp_path):
    base = signature()
    assert signature() == base
    assert ad.signature_differences(base, signature(anchor_num=5)) == ["anchor_num"]
    assert ad.signature_differences(base, signature(size=25)) == ["size"]
    assert ad.signature_differences(base, signature(random_seed=0)) == ["random_seed"]
    other = tmp_path / "other.sdf"
    other.write_text(open(LIGAND).read() + "\n")
    assert ad.signature_differences(base, signature(input_file=str(other))) == ["ligands_sha256"]
    copy = tmp_path / "copy.sdf"
    copy.write_text(open(LIGAND).read())
    assert signature(input_file=str(copy)) == base  # same content, different path: still the same run


def make_previous_run(tmp_path, **changes):
    out = tmp_path / "out"
    out.mkdir()
    (out / ad.SETTINGS_FILE).write_text(json.dumps(signature(**changes)))
    return str(out)


def test_resume_needs_an_existing_folder(tmp_path):
    with pytest.raises(ValueError, match="Nothing to resume"):
        ad.prepare_docking_directory(str(tmp_path / "missing"), False, True, signature())


def test_resume_refuses_a_folder_that_is_not_a_previous_run(tmp_path):
    (tmp_path / "stuff").mkdir()
    with pytest.raises(ValueError, match="not a previous run"):
        ad.prepare_docking_directory(str(tmp_path / "stuff"), False, True, signature())


def test_resume_refuses_changed_settings(tmp_path):
    out = make_previous_run(tmp_path, anchor_num=5)
    with pytest.raises(ValueError, match="anchor_num"):
        ad.prepare_docking_directory(out, False, True, signature())


def test_resume_accepts_identical_settings_and_keeps_the_results(tmp_path):
    out = make_previous_run(tmp_path)
    (tmp_path / "out" / "0").mkdir()
    ad.prepare_docking_directory(out, False, True, signature())
    assert (tmp_path / "out" / "0").is_dir()


def test_an_existing_folder_is_still_protected_without_resume(tmp_path):
    out = make_previous_run(tmp_path)
    with pytest.raises(ValueError, match="already exists"):
        ad.prepare_docking_directory(out, False, False)
    ad.prepare_docking_directory(out, True, False)  # --force starts again
    assert os.listdir(out) == []


def test_resume_and_force_cannot_be_combined(monkeypatch, tmp_path):
    argv = ["anchored_docking.py", "--template_file", TEMPLATE, "--input_file", LIGAND, "--receptor_file", RECEPTOR,
            "--docking_dir", str(tmp_path), "--resume", "--force"]
    monkeypatch.setattr(ad.sys, "argv", argv)
    with pytest.raises(SystemExit) as info:
        ad.parse_arguments()
    assert info.value.code == 2


# --------------------------------------------------------------------------
# Summary file
# --------------------------------------------------------------------------


def test_summary_round_trip_and_sanitising(tmp_path):
    rows = {1: {"index": 1, "name": "b", "status": "failed", "note": "line one\nline\ttwo"},
            0: {"index": 0, "name": "a", "status": "docked", "anchor_max_dev": 0.5, "core_rmsd": 1.2}}
    ad.write_summary(str(tmp_path), rows)
    read = ad.read_summary(str(tmp_path))
    assert [r["index"] for r in read] == ["0", "1"]  # ordered by ligand
    assert read[0]["anchor_max_dev"] == "0.5" and read[1]["note"] == "line one line two"
    assert not (tmp_path / "summary.tsv.tmp").exists()
    assert ad.read_summary(str(tmp_path / "nowhere")) == []


def test_failure_descriptions():
    error = subprocess.CalledProcessError(1, ["/x/mk_prepare_ligand.py", "-i", "a"])
    assert ad.describe_failure(error) == "mk_prepare_ligand.py failed (exit code 1)"
    assert ad.describe_failure(ad.LigandError("could not read")) == "could not read"
    assert ad.describe_failure(KeyError(3)).startswith("KeyError")
    assert ad.describe_failure(ValueError("a\nb")) == "ValueError: a"


# --------------------------------------------------------------------------
# The batch loop
# --------------------------------------------------------------------------


@pytest.fixture
def batch(tmp_path, monkeypatch):
    """Three ligands where the middle one cannot be read, with the per-ligand workflow stubbed out."""
    good = Chem.SDMolSupplier(LIGAND)[0]
    good.SetProp("_Name", "first")
    unreadable = Chem.MolFromSmiles("C(C)(C)(C)(C)C", sanitize=False)
    unreadable.SetProp("_Name", "second")
    last = Chem.Mol(good)
    last.SetProp("_Name", "third")
    ligands = tmp_path / "ligands.sdf"
    writer = Chem.SDWriter(str(ligands))
    for mol in (good, unreadable, last):
        writer.write(mol)
    writer.close()

    calls = []
    behaviour = {"checked": {"ok": True, "problems": [], "anchor_max": 0.4, "core_rmsd": 0.9}}

    def fake_dock(args, tools, template_mol, input_mol, docking_dir):
        calls.append(input_mol.GetProp("_Name"))
        if input_mol.GetProp("_Name") in behaviour.get("fail", ()):
            raise subprocess.CalledProcessError(2, ["autodock_gpu_128wi"])  # fails before producing a pose
        for name in ("ligand_docked.sdf", "ligand.dlg", "rec.A.map"):
            with open(os.path.join(docking_dir, name), "w") as fh:
                fh.write("x")
        with open(os.path.join(docking_dir, "pose_check.json"), "w") as fh:
            json.dump(behaviour["checked"], fh)
        return behaviour["checked"]

    monkeypatch.setattr(ad, "dock_ligand", fake_dock)
    monkeypatch.setattr(ad, "find_autodock_gpu", lambda path=None: "adgpu")
    monkeypatch.setattr(ad, "find_autogrid", lambda path=None: "autogrid4")
    out = tmp_path / "out"

    def run(*extra):
        argv = ["anchored_docking.py", "--template_file", TEMPLATE, "--input_file", str(ligands),
                "--receptor_file", RECEPTOR, "--docking_dir", str(out), *extra]
        monkeypatch.setattr(ad.sys, "argv", argv)
        ad.main()

    return run, out, calls, behaviour


def summary_by_index(out):
    return {int(r["index"]): r for r in ad.read_summary(str(out))}


def test_one_bad_ligand_does_not_stop_the_others(batch, capsys):
    run, out, calls, _ = batch
    with pytest.raises(SystemExit) as info:
        run()
    assert info.value.code == 1  # the job is marked failed because a ligand failed
    assert calls == ["first", "third"]  # the third was docked after the second failed
    rows = summary_by_index(out)
    assert [rows[i]["status"] for i in range(3)] == ["docked", "failed", "docked"]
    assert "could not read" in rows[1]["note"]
    assert (out / "1" / "error.txt").read_text().startswith("RDKit could not read")
    assert (out / "0" / "ligand_docked.sdf").is_file() and (out / "2" / "ligand_docked.sdf").is_file()
    text = capsys.readouterr().out
    assert "2 of 3 ligand(s) docked, 1 failed" in text and "--resume" in text


def test_a_crash_in_an_external_program_is_contained(batch):
    run, out, calls, behaviour = batch
    behaviour["fail"] = {"first"}
    with pytest.raises(SystemExit):
        run()
    rows = summary_by_index(out)
    assert rows[0]["status"] == "failed" and rows[0]["note"] == "autodock_gpu_128wi failed (exit code 2)"
    assert rows[2]["status"] == "docked"


def test_clean_keeps_only_results_and_the_pose_check(batch):
    run, out, _, _ = batch
    with pytest.raises(SystemExit):
        run("--clean")
    assert sorted(os.listdir(out / "0")) == ["ligand_docked.sdf", "pose_check.json"]
    assert (out / "1" / "error.txt").is_file()  # a failed ligand keeps its files for diagnosis


def test_flagged_poses_are_reported_but_not_failures(batch, capsys):
    run, out, _, behaviour = batch
    behaviour["checked"] = {"ok": False, "problems": ["an anchor atom is 4.0 A from its target"],
                            "anchor_max": 4.0, "core_rmsd": 5.1}
    with pytest.raises(SystemExit):
        run()
    rows = summary_by_index(out)
    assert rows[0]["status"] == "flagged" and rows[0]["anchor_max_dev"] == "4.0"
    assert "2 pose(s) flagged" in capsys.readouterr().out


def test_resume_skips_finished_ligands_and_retries_failed_ones(batch, capsys):
    run, out, calls, _ = batch
    with pytest.raises(SystemExit):
        run()
    calls.clear()
    with pytest.raises(SystemExit):
        run("--resume")  # the unreadable ligand fails again, so the exit code is still 1
    assert calls == []  # nothing that was already docked is run again
    assert "already docked, skipped" in capsys.readouterr().out
    rows = summary_by_index(out)
    assert [rows[i]["status"] for i in range(3)] == ["docked", "failed", "docked"]


def test_resume_redoes_a_ligand_whose_earlier_attempt_failed(batch):
    run, out, calls, behaviour = batch
    behaviour["fail"] = {"third"}
    with pytest.raises(SystemExit):
        run()
    assert (out / "2" / "error.txt").is_file()
    behaviour["fail"] = ()
    calls.clear()
    with pytest.raises(SystemExit):  # ligand 1 is still unreadable
        run("--resume")
    assert calls == ["third"]  # only the failed one was redone
    assert not (out / "2" / "error.txt").exists()
    assert summary_by_index(out)[2]["status"] == "docked"


def test_resume_does_not_mix_in_changed_settings(batch):
    run, out, _, _ = batch
    with pytest.raises(SystemExit):
        run()
    with pytest.raises(ValueError, match="anchor_num"):
        run("--resume", "--anchor_num", "5")


def test_all_good_batch_exits_cleanly(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(ad, "dock_ligand", lambda *a, **k: {"ok": True, "problems": [], "anchor_max": 0.3,
                                                          "core_rmsd": 0.5})
    monkeypatch.setattr(ad, "find_autodock_gpu", lambda path=None: "adgpu")
    monkeypatch.setattr(ad, "find_autogrid", lambda path=None: "autogrid4")
    out = tmp_path / "out"
    monkeypatch.setattr(ad.sys, "argv", ["anchored_docking.py", "--template_file", TEMPLATE, "--input_file", LIGAND,
                                         "--receptor_file", RECEPTOR, "--docking_dir", str(out)])
    ad.main()  # must not raise SystemExit
    assert "1 of 1 ligand(s) docked, 0 failed" in capsys.readouterr().out
    assert summary_by_index(out)[0]["status"] == "docked"


# --------------------------------------------------------------------------
# Spread anchor selection
# --------------------------------------------------------------------------


def test_spread_selection_picks_the_corners_of_a_square():
    square = [[0, 0, 0], [10, 0, 0], [0, 10, 0], [10, 10, 0], [5, 5, 0]]
    assert ad.spread_selection(square, 4) == [0, 1, 2, 3]
    assert ad.spread_selection(square, 2) in ([0, 3], [1, 2])      # opposite corners


def test_spread_selection_takes_the_ends_of_a_line():
    line = [[x, 0, 0] for x in range(10)]
    assert ad.spread_selection(line, 2) == [0, 9]
    assert ad.spread_selection(line, 3) == [0, 4, 9] or ad.spread_selection(line, 3) == [0, 5, 9]


def test_spread_selection_returns_everything_when_asked_for_more():
    points = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    assert ad.spread_selection(points, 17) == [0, 1, 2]
    assert ad.spread_selection(points, 0) == [] and ad.spread_selection([], 3) == []


def test_spread_selection_never_repeats_a_position():
    same = [[1.0, 1.0, 1.0]] * 4 + [[5.0, 0.0, 0.0]]
    chosen = ad.spread_selection(same, 4)
    assert len(chosen) == len(set(chosen)) == 4


def test_spread_selection_is_deterministic_and_spreads_wider_than_random_picks():
    import numpy as np

    template = Chem.SDMolSupplier(TEMPLATE)[0]
    points = template.GetConformer().GetPositions()

    def min_gap(indices):
        chosen = points[list(indices)]
        return min(np.linalg.norm(a - b) for i, a in enumerate(chosen) for b in chosen[i + 1:])

    spread = ad.spread_selection(points, 6)
    assert spread == ad.spread_selection(points, 6)
    rng = np.random.default_rng(0)
    random_gaps = [min_gap(rng.choice(len(points), 6, replace=False)) for _ in range(200)]
    assert min_gap(spread) > np.percentile(random_gaps, 95)


def anchor_search(mode, number, seed=None):
    template = Chem.SDMolSupplier(TEMPLATE)[0]
    ligand = Chem.SDMolSupplier(LIGAND)[0]
    return template, ad.find_anchor_points(template, ligand, number, 100, seed, anchor_mode=mode)


def test_spread_mode_anchors_are_core_atoms_at_their_template_positions():
    import numpy as np

    template, search = anchor_search("spread", 8)
    assert len(search.anchor_patterns) == len(search.anchor_positions) == 8
    assert len({idx for _, idx in search.anchor_patterns}) == 8
    template_positions = template.GetConformer().GetPositions()
    for _, idx in search.anchor_patterns:
        target = template_positions[search.ptn2template[idx]]
        assert any(np.allclose(target, p) for p in search.anchor_positions)


def test_spread_selection_starts_from_the_position_farthest_from_the_centre():
    # The middle point comes first in the list but must not be used when two points are asked for.
    assert ad.spread_selection([[5, 0, 0], [0, 0, 0], [10, 0, 0]], 2) == [1, 2]


def test_spread_mode_applies_the_selection_to_the_shared_core():
    import numpy as np

    _, everything = anchor_search("random", 0)
    wanted = ad.spread_selection(everything.anchor_positions, 8)
    _, search = anchor_search("spread", 8)
    assert search.anchor_patterns == [everything.anchor_patterns[i] for i in wanted]
    assert np.allclose(search.anchor_positions, [everything.anchor_positions[i] for i in wanted])


def test_spread_mode_does_not_depend_on_the_seed():
    _, first = anchor_search("spread", 8, seed=1)
    _, second = anchor_search("spread", 8, seed=2)
    assert first.anchor_patterns == second.anchor_patterns


def test_spread_mode_uses_every_usable_atom_when_the_core_is_smaller_than_the_request():
    _, everything = anchor_search("random", 0)
    available = len(everything.anchor_patterns)
    _, search = anchor_search("spread", available + 10)          # random mode would raise here
    assert len(search.anchor_patterns) == available
    with pytest.raises(ValueError, match="Not enough anchor points"):
        anchor_search("random", available + 10)


def test_spread_is_a_valid_mode_and_part_of_the_resume_signature():
    assert ad.signature_differences(signature(), signature(anchor_mode="spread")) == ["anchor_mode"]


# --------------------------------------------------------------------------
# AutoDock-GPU's atom-type limit
# --------------------------------------------------------------------------


def fake_binary(folder, limit=None, name="autodock_gpu_128wi"):
    folder.mkdir(parents=True, exist_ok=True)
    binary = folder / name
    binary.write_text("#!/bin/sh\n")
    binary.chmod(0o755)
    if limit is not None:
        (folder / ad.ATOM_TYPE_LIMIT_FILE).write_text(f"{limit}\n")
    return str(binary)


def test_a_stock_binary_has_the_stock_atom_type_limit(tmp_path):
    assert ad.atom_type_limit(fake_binary(tmp_path / "stock")) == ad.DEFAULT_ATOM_TYPE_LIMIT == 32


def test_a_custom_build_states_its_limit_beside_the_binary(tmp_path):
    assert ad.atom_type_limit(fake_binary(tmp_path / "custom", 64)) == 64


def test_the_limit_is_found_through_a_symlink(tmp_path):
    real = fake_binary(tmp_path / "custom", 48)
    link = tmp_path / "bin" / "adgpu"
    link.parent.mkdir()
    link.symlink_to(real)
    assert ad.atom_type_limit(str(link)) == 48


@pytest.mark.parametrize("content", ["not a number", "", "16", "-5"])
def test_an_unusable_limit_file_falls_back_to_the_stock_limit(tmp_path, content):
    binary = fake_binary(tmp_path / "odd")
    (tmp_path / "odd" / ad.ATOM_TYPE_LIMIT_FILE).write_text(content)
    assert ad.atom_type_limit(binary) == 32


def test_ligand_atom_types_are_read_from_the_pdbqt(tmp_path):
    pdbqt = tmp_path / "ligand.pdbqt"
    pdbqt.write_text(
        "ROOT\n"
        "ATOM      1  C01 UNL     1       0.000   0.000   0.000  0.00  0.00    +0.000 A00\n"
        "ATOM      2  C02 UNL     1       1.000   0.000   0.000  0.00  0.00    +0.000 C \n"
        "HETATM    3  O03 UNL     1       2.000   0.000   0.000  0.00  0.00    -0.300 OA\n"
        "ATOM      4  H04 UNL     1       3.000   0.000   0.000  0.00  0.00    +0.200 HD\n"
        "ATOM      5  C05 UNL     1       4.000   0.000   0.000  0.00  0.00    +0.000 C \n"
        "ENDROOT\nTORSDOF 0\n")
    assert ad.ligand_atom_types(str(pdbqt)) == {"A00", "C", "OA", "HD"}


def script_that(tmp_path, body):
    path = tmp_path / "fake_adgpu.sh"
    path.write_text("#!/bin/sh\n" + body)
    path.chmod(0o755)
    return str(path)


def test_an_error_printed_by_autodock_gpu_is_reported_even_with_exit_code_zero(tmp_path):
    """AutoDock-GPU stops with 'Too many ligand atom types' but still exits with 0; without this the run
    only failed later, in mk_export.py, hiding the cause."""
    binary = script_that(tmp_path, "echo 'Error: Too many ligand atom types (more than 32).'\n")
    with pytest.raises(ad.LigandError, match="Too many ligand atom types"):
        ad.run_docking(str(tmp_path), binary, [], {})


def test_a_missing_result_file_is_reported(tmp_path):
    with pytest.raises(ad.LigandError, match="no output file"):
        ad.run_docking(str(tmp_path), script_that(tmp_path, "echo working\n"), [], {})


def test_a_successful_run_passes_and_its_output_is_shown(tmp_path, capsys):
    binary = script_that(tmp_path, "echo progress line\ntouch ligand.dlg\n")
    ad.run_docking(str(tmp_path), binary, [], {})
    assert "progress line" in capsys.readouterr().out


def test_a_failing_exit_code_is_still_an_error(tmp_path):
    with pytest.raises(subprocess.CalledProcessError):
        ad.run_docking(str(tmp_path), script_that(tmp_path, "exit 3\n"), [], {})


needs_meeko = pytest.mark.skipif(shutil.which("mk_prepare_ligand.py") is None,
                                 reason="Meeko is not on PATH (run the tests through pixi)")


def prepare_with_limit(tmp_path, limit, mode="spread", number=22, indices=None):
    """The real anchor choice and Meeko preparation for the crystal ligand docked onto itself (22 usable anchors)."""
    import argparse

    template = Chem.SDMolSupplier(TEMPLATE)[0]
    args = argparse.Namespace(anchor_num=number, samples=100, random_seed=None, query_smarts=None,
                              anchor_indices=indices, anchor_mode=mode, verbose=False)
    tools = {"scripts_path": os.path.join(REPO, "scripts"), "atom_type_limit": limit}
    result = ad.prepare_anchored_ligand(args, tools, template, template, str(tmp_path))
    return result, ad.ligand_atom_types(str(tmp_path / "ligand.pdbqt"))


@needs_meeko
def test_all_anchors_are_used_when_the_ligand_fits(tmp_path):
    (search, _, _, note), types = prepare_with_limit(tmp_path, 32)
    assert len(search.anchor_patterns) == 22 and note == "" and len(types) <= 32


@needs_meeko
def test_spread_mode_drops_anchors_until_the_ligand_fits(tmp_path):
    (search, anchors, mapping, note), types = prepare_with_limit(tmp_path, 20)
    assert len(types) <= 20
    assert 1 <= len(search.anchor_patterns) < 22
    assert len(anchors) == len(mapping) == len(search.anchor_patterns)
    assert f"only {len(search.anchor_patterns)} of 22 anchors" in note and "20 atom types" in note
    # the files on disk are those of the final choice, not of the first attempt
    assert len(json.load(open(tmp_path / "anchor_positions.json"))) == len(search.anchor_patterns)
    base = json.load(open(os.path.join(REPO, "scripts", "base_parameters.json")))["ad4"]
    assert len(json.load(open(tmp_path / "parameters.json"))["ad4"]) - len(base) == len(search.anchor_patterns)


@needs_meeko
def test_a_ligand_that_needs_exactly_the_limit_keeps_all_its_anchors(tmp_path):
    for folder in ("probe", "exact", "less"):
        (tmp_path / folder).mkdir()
    _, types = prepare_with_limit(tmp_path / "probe", 100)        # all 22 anchors: how many types do they need?
    needed = len(types)
    (search, _, _, note), _ = prepare_with_limit(tmp_path / "exact", needed)
    assert len(search.anchor_patterns) == 22 and note == ""
    (reduced, _, _, note), _ = prepare_with_limit(tmp_path / "less", needed - 1)
    assert len(reduced.anchor_patterns) < 22 and note


@needs_meeko
def test_the_number_of_anchors_is_lowered_by_the_excess_not_one_at_a_time(tmp_path, monkeypatch):
    calls = []
    real = ad.prepare_ligand
    monkeypatch.setattr(ad, "prepare_ligand", lambda *a, **k: (calls.append(1), real(*a, **k))[1])
    prepare_with_limit(tmp_path, 14)                 # 22 anchors asked for, far more than fit
    assert len(calls) <= 6                            # measured: 4; one at a time would take more than a dozen


@needs_meeko
@pytest.mark.parametrize("mode, number, limit", [("random", 22, 14), ("area", 3, 4)])
def test_other_modes_fail_with_an_explanation_instead(tmp_path, mode, number, limit):
    with pytest.raises(ad.LigandError, match="atom types.*spread mode"):
        prepare_with_limit(tmp_path, limit, mode=mode, number=number)


@needs_meeko
def test_a_ligand_that_cannot_fit_even_one_anchor_fails(tmp_path):
    with pytest.raises(ad.LigandError, match="atom types"):
        prepare_with_limit(tmp_path, 2)


@needs_meeko
def test_a_higher_limit_allows_the_full_set_where_the_stock_limit_does_not(tmp_path):
    big = Chem.MolFromSmiles("CC1=C(C(=CC=C1)C)OCC(=O)NC(CC2=CC=CC=C2)C(CC(CC3=CC=CC=C3)NC(=O)C(C(C)C)N4CCCNC4=O)O")
    big = Chem.AddHs(big)
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(big, randomSeed=3)
    import argparse
    args = argparse.Namespace(anchor_num=64, samples=100, random_seed=None, query_smarts=None,
                              anchor_indices=None, anchor_mode="spread", verbose=False)
    scripts = os.path.join(REPO, "scripts")
    stock, custom = tmp_path / "stock", tmp_path / "custom"
    stock.mkdir()
    custom.mkdir()
    _, _, _, stock_note = ad.prepare_anchored_ligand(args, {"scripts_path": scripts, "atom_type_limit": 32}, big, big, str(stock))
    search, _, _, custom_note = ad.prepare_anchored_ligand(args, {"scripts_path": scripts, "atom_type_limit": 64}, big, big, str(custom))
    assert "only" in stock_note and custom_note == ""             # 30 usable anchors do not fit 32 types, but fit 64
    assert len(search.anchor_patterns) == 30
    assert len(ad.ligand_atom_types(str(stock / "ligand.pdbqt"))) <= 32


def test_a_note_about_dropped_anchors_reaches_the_summary_and_the_log_line():
    checked = {"ok": True, "problems": [], "anchor_max": 0.9, "core_rmsd": 0.6,
               "limit_note": "only 20 of 22 anchors could be used"}
    row = ad.row_from_check(0, "x", checked)
    assert row["status"] == "docked" and row["note"] == "only 20 of 22 anchors could be used"
    assert "anchor deviation 0.9 A" in ad.describe_row(row) and "only 20 of 22 anchors" in ad.describe_row(row)
    plain = ad.row_from_check(0, "x", {"ok": True, "problems": [], "anchor_max": 0.9, "core_rmsd": 0.6})
    assert plain["note"] == "" and ad.describe_row(plain).endswith("core RMSD 0.6 A)")
