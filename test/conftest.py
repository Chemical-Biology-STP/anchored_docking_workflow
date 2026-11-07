import os
import sys

import pytest

REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_DIR not in sys.path:
    sys.path.insert(0, REPO_DIR)


def propka_block(pdb_text):
    """A PROPKA summary like the one in a real PDB2PQR log, with a made-up pKa for each titratable residue."""
    typical = {"ASP": 3.8, "GLU": 4.5, "LYS": 10.5, "ARG": 12.5, "HIS": 6.5, "CYS": 9.0, "TYR": 10.0}
    seen, rows = set(), []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            key = (line[17:20], line[21:22], int(line[22:26]))
            if key not in seen and key[0] in typical:
                seen.add(key)
                rows.append(f"   {key[0]} {key[2]:>3} {key[1]}  {typical[key[0]]:>7.2f}      {typical[key[0]]:.2f}")
    return "SUMMARY OF THIS PREDICTION\n       Group      pKa  model-pKa   ligand atom-type\n" + "\n".join(rows) + "\n-----\n"


class FakePdb2pqr:
    """Stands in for PDB2PQR: adds one hydrogen to the first residue and writes a PROPKA-style log."""

    def __init__(self):
        self.calls = []
        self.fail = False

    def __call__(self, pdb_text, ph):
        import protein_prep as pp

        self.calls.append(ph)
        if self.fail:
            raise pp.ProteinPrepError("PDB2PQR failed (exit code 1). pretend error")
        out = []
        added = False
        for line in pdb_text.splitlines():
            out.append(line)
            if line.startswith("ATOM") and not added:
                padded = line.ljust(80)
                out.append(padded[:12] + " H1 " + padded[16:76] + " H" + padded[78:])
                added = True
        return "\n".join(out) + "\n", propka_block(pdb_text)


@pytest.fixture
def fake_pdb2pqr(monkeypatch):
    import protein_prep as pp

    fake = FakePdb2pqr()
    monkeypatch.setattr(pp, "run_pdb2pqr", fake)
    return fake


@pytest.fixture(autouse=True)
def scheduler_allows_everything(monkeypatch):
    """The wizard asks Slurm whether the account may use a partition; tests must not depend on who runs them.
    Tests of that check itself call dock_cli.check_submission through the original function they imported."""
    import dock_cli

    monkeypatch.setattr(dock_cli, "check_submission", lambda *args, **kwargs: dock_cli.SubmitCheck(True))
