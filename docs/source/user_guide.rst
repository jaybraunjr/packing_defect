User Guide
==========

This guide walks through common workflows for analyzing membrane packing defects with ``packing_defect``.

Prerequisites
-------------

- MDAnalysis installed (use the conda env in ``devtools/conda-envs/test_env.yaml``)
- Topology and trajectory that MDAnalysis can read (e.g., GRO/XTC, PSF/DCD)


Typical Workflow: Packing Analysis
----------------------------------

Use the topology + trajectory path when your simulation contains lipids of interest and you want the toolkit to classify atoms and stamp radii automatically::

   python -m packing_defect.run_defect \
     --top system.gro \
     --traj traj.xtc \
     --out results \
     --leaflet both \
     --start 0 --stop 100 --stride 5

Key options
- ``--leaflet {both,up,dw}``: restricts analysis to selected leaflet(s)
- ``--start/--stop/--stride``: frame slicing over the trajectory
- ``--class-json``: optional path to a custom classifier JSON to override defaults

Outputs
- Per-defect-type subfolders under ``results/`` containing ``*_frame_<i>.gro``
- ``<TYPE>.dat`` two-column size distributions (bin center, probability)
- Optional plots via the analyzer's ``plot()`` method


Radius Workflow: Pre-split GRO Frames
-------------------------------------

When you already have GRO frames split per lipid type and want to filter by protein proximity and measure defect sizes::

   python -m packing_defect.run_radius \
     --input inputs \
     --output results_radius \
     --lipids PLacyl TGacyl TGglyc \
     --start 0 --end 100 \
     --protein-count 627 \
     --cutoff 1.5

Notes
- ``--protein-count`` is the number of protein atoms at the top of each GRO
- ``--no-cutoff`` disables the protein-distance filter


Interpreting Results
--------------------

- ``*.dat`` files are normalized histograms of connected-component sizes on a 2D grid
- Sizes are in units of grid cells; choose grid spacing via analyzer parameters (default 1.0 Å)
- Per-leaflet results are aggregated into "up"/"dw" and "combined" where applicable


Extensibility
-------------

- High-level analyzers live in :mod:`packing_defect.core.analyzers`
- The 2D grid accumulator is :class:`packing_defect.core.grid.DefectGrid`
- Clustering utilities are in :mod:`packing_defect.core.cluster`

Refer to the :doc:`api` for full reference documentation.

