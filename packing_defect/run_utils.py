"""Small helper(s) to run analyzers consistently from CLIs.

This module centralizes the typical sequence used by the run_* entrypoints
so different analyzers behave the same way when invoked as scripts.
"""

from typing import Protocol


class _RunnableAnalyzer(Protocol):
    def run(self) -> None: ...
    def plot(self) -> None: ...
    def save_results(self, filename: str, data) -> str: ...
    results: dict


def run_analysis(analyzer: _RunnableAnalyzer) -> None:
    """Execute, plot, and persist results for a defect analyzer.

    The standard sequence used by CLI entrypoints:
    1) ``run()``, 2) ``plot()``, 3) ``save_results()`` to ``results.txt``.

    Parameters
    ----------
    analyzer : BaseDefectAnalyzer-like
        Any object exposing ``run()``, ``plot()``, ``save_results()`` and
        a ``results`` attribute.
    """
    analyzer.run()
    analyzer.plot()
    analyzer.save_results("results.txt", analyzer.results)
