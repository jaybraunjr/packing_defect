from pathlib import Path
from typing import Dict, Tuple, Callable, Union

# If you have a ClassificationStrategy class, import it for isinstance checks
from packing_defect.core.classification import ClassificationStrategy


class TopologyReader:
    """
    Read CHARMM RTF/STR topology files (``*.rtf``, ``*.str``) and extract
    atom properties for a given residue.

    Each atom is mapped to a tuple ``(radius, code)``, where:

    - **radius** is looked up from a provided radii dictionary
    - **code** is determined using a classification strategy

    The classifier may be either:

    * An object with a ``.classify(resname, atom_name)`` method
    * A plain function ``f(resname, atom_name) -> int``
    """

    def __init__(
        self,
        radii: Dict[str, float],
        classifier: Union[ClassificationStrategy, Callable[[str, str], int]],
    ):
        """
        Parameters
        ----------
        radii : dict of {str: float}
            Mapping from atom type to radius.
        classifier : ClassificationStrategy or callable
            Either a strategy instance with ``classify(resname, atom_name)``
            or a simple function ``(resname, atom_name) -> int`` that
            returns a classification code.
        """
        self.radii = radii
        self.classifier = classifier

    def read(self, resname: str, topo_path: str) -> Dict[str, Tuple[float, int]]:
        """
        Parse the ``RESI`` block for a given residue in the topology file.

        Parameters
        ----------
        resname : str
            Residue name to search for (e.g., ``POPC``).
        topo_path : str
            Path to a CHARMM topology file (``.rtf`` or ``.str``).

        Returns
        -------
        dict of {str: (float, int)}
            Mapping from atom name to a tuple ``(radius, code)``.

        Raises
        ------
        FileNotFoundError
            If the topology file does not exist.
        KeyError
            If an atom type is not found in the radii map.
        TypeError
            If the classifier is not callable or does not implement
            ``.classify()``.
        """
        path = Path(topo_path)
        if not path.exists():
            raise FileNotFoundError(f"Topology file not found: {topo_path}")

        output: Dict[str, Tuple[float, int]] = {}
        in_resi = False

        with path.open() as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("!"):
                    continue

                # Start when we see "RESI <resname>"
                if line.upper().startswith(f"RESI {resname.upper()}"):
                    in_resi = True
                    continue

                # End when BOND block starts
                if in_resi and line.upper().startswith("BOND"):
                    break

                # Parse ATOM lines
                if in_resi and line.upper().startswith("ATOM"):
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    atom_name, atom_type = parts[1], parts[2]

                    # Look up radius
                    if atom_type not in self.radii:
                        raise KeyError(f"Atom type '{atom_type}' not in radii map")
                    radius = float(self.radii[atom_type])

                    # Classify: support both objects and bare functions
                    if hasattr(self.classifier, "classify"):
                        code = self.classifier.classify(resname, atom_name)
                    elif callable(self.classifier):
                        code = self.classifier(resname, atom_name)
                    else:
                        raise TypeError("Classifier must be a function or have a .classify()")

                    output[atom_name] = (radius, code)

        return output
