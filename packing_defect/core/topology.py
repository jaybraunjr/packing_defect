# core/topology.py

from pathlib import Path
from typing import Dict, Tuple, Callable, Union

# If you have a ClassificationStrategy class, import it for isinstance checks
from packing_defect.core.classification import ClassificationStrategy

class TopologyReader:
    """
    Reads a CHARMM RTF/STR topology (*.rtf, *.str) and for each atom
    returns (radius, code). The classifier can be either:
      - An object with a .classify(resname, atom_name) method, OR
      - A plain function f(resname, atom_name) -> int
    """
    def __init__(self,
                 radii: Dict[str, float],
                 classifier: Union[ClassificationStrategy, Callable[[str,str],int]]):
        """
        Parameters
        ----------
        radii : dict mapping atom_type -> radius (float)
        classifier : either a ClassificationStrategy instance, or a function
                     taking (resname, atom_name) and returning an int code
        """
        self.radii = radii
        self.classifier = classifier

    def read(self, resname: str, topo_path: str) -> Dict[str, Tuple[float,int]]:
        """
        Parse the RESI block for `resname` in the given topology file.
        Returns { atom_name: (radius, code) }.
        """
        path = Path(topo_path)
        if not path.exists():
            raise FileNotFoundError(f"Topology file not found: {topo_path}")

        output: Dict[str, Tuple[float,int]] = {}
        in_resi = False

        with path.open() as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('!'):
                    continue

                # start when we see RESI <resname>
                if line.upper().startswith(f'RESI {resname.upper()}'):
                    in_resi = True
                    continue

                # end when BOND block starts
                if in_resi and line.upper().startswith('BOND'):
                    break

                # parse ATOM lines
                if in_resi and line.upper().startswith('ATOM'):
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    atom_name, atom_type = parts[1], parts[2]

                    # look up radius
                    if atom_type not in self.radii:
                        raise KeyError(f"Atom type '{atom_type}' not in radii map")
                    radius = float(self.radii[atom_type])

                    # classify: support both objects and bare functions
                    if hasattr(self.classifier, 'classify'):
                        code = self.classifier.classify(resname, atom_name)
                    elif callable(self.classifier):
                        code = self.classifier(resname, atom_name)
                    else:
                        raise TypeError("Classifier must be a function or have a .classify()")

                    output[atom_name] = (radius, code)

        return output
