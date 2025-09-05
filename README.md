<img src="https://github.com/user-attachments/assets/32cbe599-252c-4f6d-b1e1-a49f109eb614" width="750" />

[//]: # (Badges)

| **Docs** | [![Documentation Status][badge_docs]][url_docs] |
| :-- | :-- |
| **CI** | [![GH Actions Status][badge_actions]][url_actions] [![codecov][badge_codecov]][url_codecov] |
| **Project** | [![License: GPL v2][badge_license]][url_license]  [![Powered by MDAnalysis][badge_mda]][url_mda] |

[badge_actions]: https://github.com/jaybraunjr/packing_defect/actions/workflows/gh-ci.yaml/badge.svg
[badge_codecov]: https://codecov.io/gh/jaybraunjr/packing_defect/branch/main/graph/badge.svg
[badge_docs]: https://readthedocs.org/projects/packing-defect/badge/?version=latest
[badge_license]: https://img.shields.io/badge/License-GPLv2-blue.svg
[badge_mda]: https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg
[url_actions]: https://github.com/jaybraunjr/packing_defect/actions?query=workflow%3Agh-ci
[url_codecov]: https://codecov.io/gh/jaybraunjr/packing_defect/branch/main
[url_docs]: https://packing-defect.readthedocs.io/en/latest/
[url_license]: https://www.gnu.org/licenses/gpl-2.0
[url_mda]: https://www.mdanalysis.org


Packing Defect Analysis
-----------------------

Analyze membrane packing defects from molecular dynamics trajectories using grid‑based stamping and connected‑component clustering. The toolkit classifies atoms into defect types, writes per‑frame outputs, and summarizes defect size distributions.

Highlights
- Detect defects on upper/lower leaflets per frame
- Classify atoms by type (e.g., PL tails, TG glycerol, TG tails)
- Export per‑class `.gro` frames and `.dat` distributions
- Plot quick defect summaries and histograms
- Optional radius‑based filtering and renumbering for `.gro` frames


Install
- Conda (recommended):
  - `mamba env update --name packing_defect --file devtools/conda-envs/test_env.yaml`
  - `python -m pip install .`
- Pip:
  - `pip install .` or `pip install ".[test,doc]"` for dev extras


Quickstart
- Packing analysis (topology + trajectory):
  - `python -m packing_defect.run_defect --top system.gro --traj traj.xtc --out results --leaflet both`
- Radius workflow over pre‑split GRO frames:
  - `python -m packing_defect.run_radius --input inputs --output results_radius --lipids PLacyl TGacyl TGglyc --start 0 --end 100 --protein-count 627 --cutoff 1.5`

Full guide and API: [packing‑defect.readthedocs.io][url_docs]


License and Citation
- GPLv2; see `LICENSE`.
- Please cite [MDAnalysis](https://github.com/MDAnalysis/mdanalysis#citation) when using this toolkit.
