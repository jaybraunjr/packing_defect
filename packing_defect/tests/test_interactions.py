import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from packing_defect.core.interactions import ProteinDefectInteractions

def _make_universe_with_protein_and_others():
    n = 4
    u = mda.Universe.empty(
        n_atoms=n, n_residues=n, n_segments=1, trajectory=True,
        atom_resindex=np.arange(n), residue_segindex=np.zeros(n, int)
    )
    # Standard protein residue names so "protein" selection works
    u.add_TopologyAttr("resname", ["ALA", "ALA", "OTH", "OTH"])
    u.add_TopologyAttr("name", ["CA", "CB", "X", "Y"])
    u.add_TopologyAttr("resid", [1, 2, 3, 4])

    u.add_TopologyAttr("masses", np.ones(n))  # any positive values are fine

    coords = np.array([
        [0.0, 0.0, 0.0],   # protein (ALA)
        [0.0, 0.5, 0.0],   # protein (ALA)
        [0.1, 0.1, 0.0],   # close “defect”
        [10.0, 10.0, 0.0]  # far
    ])
    u.trajectory = MemoryReader(coords[np.newaxis, :, :])
    return u

def test_interactions_detects_hits_within_cutoff():
    u = _make_universe_with_protein_and_others()
    tracker = ProteinDefectInteractions(cutoff=0.3)
    tracker.track_frame(u, frame_index=0)
    assert 0 in tracker.results
    hits = tracker.results[0]
    assert any(resname == "ALA" for (_resid, resname) in hits)
