import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from pathlib import Path
import pytest

from packing_defect.core.topology import TopologyReader


def _write_topology(tmp_path: Path):
    txt = """
* Simple CHARMM RTF-like
RESI POPC 0.00
ATOM C1   CTL2  0.00
ATOM P    P1    0.00
ATOM C2   CTL2  0.00
BOND C1  C2
""".strip()
    p = tmp_path / "test.rtf"
    p.write_text(txt)
    return p


def test_topo_read(tmp_path: Path):
    topo = _write_topology(tmp_path)
    radii = {"CTL2": 1.2, "P1": 1.5}

    # Classifier function: return 9 for everything
    def classify(res, atom):
        return 9

    reader = TopologyReader(radii=radii, classifier=classify)
    out = reader.read("POPC", str(topo))
    assert out["C1"] == (1.2, 9)
    assert out["P"] == (1.5, 9)
    assert out["C2"] == (1.2, 9)


def test_topo_missing_file(tmp_path: Path):
    reader = TopologyReader(radii={}, classifier=lambda r, a: 0)
    with pytest.raises(FileNotFoundError):
        reader.read("POPC", str(tmp_path / "nope.rtf"))


def test_topo_missing_type(tmp_path: Path):
    topo = _write_topology(tmp_path)
    reader = TopologyReader(radii={"P1": 1.5}, classifier=lambda r, a: 1)
    with pytest.raises(KeyError):
        reader.read("POPC", str(topo))


def test_topo_bad_classifier(tmp_path: Path):
    topo = _write_topology(tmp_path)
    class NotCallable:
        pass
    bad = NotCallable()
    reader = TopologyReader(radii={"CTL2": 1.0, "P1": 1.0}, classifier=bad)
    with pytest.raises(TypeError):
        reader.read("POPC", str(topo))
