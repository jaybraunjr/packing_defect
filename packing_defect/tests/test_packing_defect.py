import json
import os
import numpy as np
from packing_defect.core.classification import DefaultClassification, UserDictClassification
from packing_defect.core.topology import TopologyReader


def test_default_classification():
    classifier = DefaultClassification()
    assert classifier.classify('POPC', 'C216') == 1
    assert classifier.classify('POPC', 'H1') == -1
    assert classifier.classify('TRIO', 'O11') == 2
    assert classifier.classify('TRIO', 'XYZ') == 3

def test_userdict_classification(tmp_path):
    json_file = tmp_path / "rules.json"
    with open(json_file, 'w') as f:
        json.dump({"RES": {"A": 9}}, f)
    classifier = UserDictClassification.from_json(str(json_file))
    assert classifier.classify('RES', 'A') == 9
    assert classifier.classify('RES', 'B') == -1


def test_topology_reader(tmp_path):
    # simple topology file with two atoms of type CTL2
    top = tmp_path / "test.top"
    top.write_text("RESI POPC 0.00\nATOM C216 CTL2\nATOM H1 CTL2\nBOND C216 H1\n")
    radii = {"CTL2": 1.5}
    reader = TopologyReader(radii, DefaultClassification())
    result = reader.read('POPC', str(top))
    assert result["C216"] == (1.5, 1)
    assert result["H1"] == (1.5, -1)
