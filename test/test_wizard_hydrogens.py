"""The wizard's handling of hydrogens for display on the 'start from a complex' route."""

import pytest
from rdkit import Chem

import structure_prep as sp
import test_results as tr
import test_structure_prep as tsp

LIGAND = tsp.LIGAND


@pytest.fixture(autouse=True)
def clear_warnings():
    tsp.d.WARNINGS.clear()
    yield
    tsp.d.WARNINGS.clear()


def template_dictionary():
    return tr.ccd_for(Chem.SDMolSupplier(tsp.TEMPLATE)[0], comp="LIG")


def run(monkeypatch, tmp_path, fetch):
    monkeypatch.setattr(sp, "fetch_ccd", fetch)
    structure = tsp.write_complex(tmp_path)
    tsp.run_complex_wizard(monkeypatch, ["1", structure, "", "no", LIGAND, "", "", *tsp.tail(tmp_path / "w" / "o")])
    return tmp_path / "w"


def test_template_hydrogens_come_from_the_dictionary_when_it_can_be_reached(monkeypatch, tmp_path, capsys):
    asked = []

    def fetch(name, timeout=15):
        asked.append(name)
        return template_dictionary()

    work = run(monkeypatch, tmp_path, fetch)
    text = capsys.readouterr().out
    assert asked.count("LIG") == 1                       # the template's entry is reused, not downloaded a second time
    assert "Template hydrogens for display: added from the component dictionary" in text
    view = work / "dock_CAT-13f_template_with_hydrogens.sdf"
    mol = Chem.SDMolSupplier(str(view), removeHs=False)[0]
    assert sum(a.GetAtomicNum() == 1 for a in mol.GetAtoms()) == 19
    assert f"--template_view_file {view}" in (work / "submit_dock_CAT-13f.sh").read_text()
    # the docking itself still uses the heavy-atom template
    docking = Chem.SDMolSupplier(str(work / "dock_CAT-13f_template.sdf"), removeHs=False)[0]
    assert sum(a.GetAtomicNum() == 1 for a in docking.GetAtoms()) == 0


def test_without_the_dictionary_the_template_is_shown_without_guessed_hydrogens(monkeypatch, tmp_path, capsys):
    work = run(monkeypatch, tmp_path, lambda name, timeout=15: None)
    text = capsys.readouterr().out
    assert "Template hydrogens for display: not added" in text and "could not be fetched" in text
    view = work / "dock_CAT-13f_template_with_hydrogens.sdf"
    mol = Chem.SDMolSupplier(str(view), removeHs=False)[0]
    assert sum(a.GetAtomicNum() == 1 for a in mol.GetAtoms()) == 0
    assert f"--template_view_file {view}" in (work / "submit_dock_CAT-13f.sh").read_text()  # still passed, so nothing is guessed


def test_a_mismatched_dictionary_entry_is_explained(monkeypatch, tmp_path, capsys):
    run(monkeypatch, tmp_path, lambda name, timeout=15: tr.ccd_for(Chem.MolFromSmiles("CCO"), comp="LIG"))
    text = capsys.readouterr().out
    assert "Template hydrogens for display: not added" in text and "does not match the component dictionary" in text


def test_the_dictionary_is_asked_once_for_each_molecule_listed_and_for_nothing_else(monkeypatch, tmp_path):
    asked = []
    run(monkeypatch, tmp_path, lambda name, timeout=15: asked.append(name))
    # the molecules in the template table (so the user can tell the drug from glycerol and sulfate), each once;
    # the zinc ion is too small to be a candidate and is never looked up
    assert sorted(asked) == ["GOL", "LIG", "SO4"]
