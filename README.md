<img src="https://github.com/user-attachments/assets/32cbe599-252c-4f6d-b1e1-a49f109eb614" width="750" />



## ** ⚠️ Under Development**

**This package is actively being developed and will soon become part of [MDAkit](https://www.mdanalysis.org/pages/mdakits/).**  
Interfaces, outputs, and command-line options are still subject to change.




[//]: # (Badges)

| **Latest release** | [![Last release tag][badge_release]][url_latest_release] ![GitHub commits since latest release (by date) for a branch][badge_commits_since]  [![Documentation Status][badge_docs]][url_docs]|
| :----------------- | :------- |
| **Status**         | [![GH Actions Status][badge_actions]][url_actions] [![codecov][badge_codecov]][url_codecov] |
| **Community**      | [![License: GPL v2][badge_license]][url_license]  [![Powered by MDAnalysis][badge_mda]][url_mda]|

[badge_actions]: https://github.com/jaybraunjr/packing_defect/actions/workflows/gh-ci.yaml/badge.svg
[badge_codecov]: https://codecov.io/gh/jaybraunjr/packing_defect/branch/main/graph/badge.svg
[badge_commits_since]: https://img.shields.io/github/commits-since/jaybraunjr/packing_defect/latest
[badge_docs]: https://readthedocs.org/projects/packing_defect/badge/?version=latest
[badge_license]: https://img.shields.io/badge/License-GPLv2-blue.svg
[badge_mda]: https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
[badge_release]: https://img.shields.io/github/release-pre/jaybraunjr/packing_defect.svg
[url_actions]: https://github.com/jaybraunjr/packing_defect/actions?query=branch%3Amain+workflow%3Agh-ci
[url_codecov]: https://codecov.io/gh/jaybraunjr/packing_defect/branch/main
[url_docs]: https://packing_defect.readthedocs.io/en/latest/?badge=latest
[url_latest_release]: https://github.com/jaybraunjr/packing_defect/releases
[url_license]: https://www.gnu.org/licenses/gpl-2.0
[url_mda]: https://www.mdanalysis.org

# Membrane packing defect analysis

This toolkit analyzes membrane packing defects from molecular dynamics simulations using grid-based detection and clustering. 
It supports defect classification by atom types and radius-based stamping, outputs frames with defects, and computes defect size distributions and statistics.

- Detect defects per frame on membrane leaflets (upper/lower)
- Classify atoms into defect types (e.g., PL tails, TG glycerol, TG tails)
- Write `.gro` or `.pdb` files for defects
- Generate `.dat` defect size distributions
- Visualize defect size histograms
- Filter/renumber `.gro` or `.pdb` files
- Compute summary statistics of large defect clusters


## Usage Examples

### 1. Run full radius-based defect analysis

```bash
python -m packing_defect.run_defect \
  --top system.gro \
  --traj trajectory.xtc \
  --out resultsPLacyl \
  --leaflet both
```

**Optional arguments**:
- `--class my_classifier.json`: use a custom classification JSON
- `--json-only`: skip topology parsing and only use the classification JSON
- `--leaflet up|dw|both`: restrict defect detection to one leaflet

This writes `.gro` files and `.dat` distributions (e.g., `PLacyl.dat`, `TGacyl.dat`, etc.) into the output directory.

---

### 2. Count large defect clusters per frame

```bash
python -m packing_defect.core.stats resultsPLacyl -c 8 --summary
```

- `-c 8` sets the minimum size for a cluster (in grid cells)
- `--summary` prints mean, median, and standard deviation across frames

---

### 3. Visualize defect size distributions

```bash
python -m packing_defect.vis resultsPLacyl \
  --title "Defect Size Distribution" \
  -o plots/PLacyl_defects.png
```

- Plots data from `PLacyl.dat`, `TGglyc.dat`, `TGacyl.dat`
- Y-axis is log-scaled for visibility of rare large defects

---

### 4. Run simplified radius stamping + filtering + histogram (optional)

```bash
python -m packing_defect.run_radius \
  --input outputs_custom4 \
  --output results_radius \
  --lipids PLacyl TGacyl TGglyc \
  --start 0 \
  --end 100 \
  --protein-count 627 \
  --cutoff 1.5
```

- Applies a distance cutoff to exclude defect atoms too close to protein
- Outputs corrected and renumbered `.gro` files
- Generates log-scale histograms (`combined_defect_histogram2.png`, etc.)

---

### 5. Track protein-defect interactions

```bash
python -m packing_defect.core.interactions \
  -i resultsPLacyl resultsTGacyl resultsTGglyc \
  -o results/interactions.csv \
  -c 1.5
```

Checks whether any protein residues come within `1.5 Å` of defect atoms in each frame and writes the interactions to a CSV file.



### Installation

#### With conda

Create a virtual environment and activate it:

```
conda create --name packing_defect
conda activate packing_defect
```
Install the development and documentation dependencies:
```
conda env update --name packing_defect --file devtools/conda-envs/test_env.yaml
conda env update --name packing_defect --file docs/requirements.yaml
```
Build this package from source:
```
pip install -e .
```
If you want to update your dependencies (which can be risky!), run:

```
conda update --all
```
And when you are finished, you can exit the virtual environment with:
```
conda deactivate
```

#### With pip
To build the package from source, run:
```
pip install .
```

If you want to create a development environment, install
the dependencies required for tests and docs with:
```
pip install ".[test,doc]"
```

### Copyright

The packing_defect source code is hosted at https://github.com/jaybraunjr/packing_defect
and is available under the GNU General Public License, version 2 (see the file [LICENSE](https://github.com/jaybraunjr/packing_defect/blob/main/LICENSE)).

Copyright (c) 2025, R. Jay Braun


#### Acknowledgements
 

Please cite [MDAnalysis](https://github.com/MDAnalysis/mdanalysis#citation) when using packing_defect in published work.
