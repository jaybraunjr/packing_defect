# packing_defect/core/analyzers/base.py

from abc import ABC, abstractmethod
import os
from dataclasses import dataclass
from typing import Any, Dict

class BaseDefectAnalyzer(ABC):
    """
    Abstract base class for all defect analyzers.
    Defines a standard interface for running analyses and plotting results.
    """

    def __init__(self, universe, output_dir: str, lipid_types=None):
        """
        Parameters
        ----------
        universe : MDAnalysis.Universe
            MDAnalysis universe containing system topology and trajectory.
        output_dir : str
            Path to directory where results should be written.
        lipid_types : list[str], optional
            List of lipid residue names to analyze.
        """
        self.universe = universe
        self.output_dir = output_dir
        self.lipid_types = lipid_types or []
        self.results = {} 
        os.makedirs(self.output_dir, exist_ok=True)

    @abstractmethod
    def run(self):
        """
        Execute the analysis over the trajectory.
        Must populate `self.results`.
        """
        pass

    @abstractmethod
    def plot(self, *args, **kwargs):
        """
        Generate plots of the analysis results.
        """
        pass


    def save_results(self, filename: str, data):
        """
        Save results to a file inside `output_dir`.
        """
        path = os.path.join(self.output_dir, filename)
        with open(path, "w") as f:
            f.write(str(data))
        print(f"[BaseDefectAnalyzer] Results written to {path}")
        return path



@dataclass
class AnalysisResult:
    data: Dict[str, Any]  # e.g., histograms, masks, coordinates
    meta: Dict[str, Any]  # e.g., params, frames processed