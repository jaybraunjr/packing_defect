Getting Started
===============

Install
-------

We recommend a conda environment for MDAnalysis and scientific deps::

   mamba env update --name packing_defect --file devtools/conda-envs/test_env.yaml
   conda activate packing_defect
   python -m pip install .

Alternatively with pip::

   pip install .
   # for tests + docs during development
   pip install ".[test,doc]"


Quickstart
----------

Run a packing analysis using a topology and trajectory::

   python -m packing_defect.run_defect \
     --top system.gro \
     --traj traj.xtc \
     --out results \
     --leaflet both

Radius workflow over pre-split GRO frames::

   python -m packing_defect.run_radius \
     --input inputs \
     --output results_radius \
     --lipids PLacyl TGacyl TGglyc \
     --start 0 \
     --end 100 \
     --protein-count 627 \
     --cutoff 1.5


Next Steps
----------

- Explore the API reference for analyzers and helpers: :doc:`api`
- See CLI entrypoints in modules: :mod:`packing_defect.run_defect`, :mod:`packing_defect.run_radius`
