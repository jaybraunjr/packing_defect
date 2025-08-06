#!/usr/bin/env python3
"""
vs.py

Plot defect size distributions for one run of packing_defect.
"""
import matplotlib
matplotlib.use('Agg')
import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


DEFAULT_LABELS  = ['PL-acyl', 'TG-glyc', 'TG-acyl']
DEFAULT_COLORS  = ['blue',    'red',      'green']
DEFAULT_MARKERS = ['o',       's',        '^']


def plot_defect_data(file_paths, labels, colors, markers, title=None, output_path=None):
    """
    Scatter defect size vs probability on a log-y axis.
    file_paths : list of three paths in order [TGacyl, TGglyc, PLacyl] or any order matching labels.
    """
    fig, ax = plt.subplots(figsize=(6,4))

    for fp, lbl, mkr in zip(file_paths, labels, markers):
        data = np.loadtxt(fp)
        ax.scatter(data[:,0], data[:,1], label=lbl, marker=mkr, s=20, alpha=0.7)

    ax.set_yscale('log')
    ax.set_ylim(bottom=1e-5)
    ax.set_xlabel('Defect size (Å²)')
    ax.set_ylabel('Probability')
    if title:
        ax.set_title(title)
    ax.legend()
    fig.tight_layout()

    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path)
        print(f'Plot saved to {output_path}')
    else:
        plt.show()

def main():
    p = argparse.ArgumentParser(description="Visualize packing_defect .dat output")
    p.add_argument('run_dir',
        help="Directory containing TGacyl.dat, TGglyc.dat, PLacyl.dat")
    p.add_argument('-o','--out',
        help="Where to save the figure (e.g. plots/vs.png). If omitted, shows interactively.")
    p.add_argument('--labels', nargs=3, default=DEFAULT_LABELS,
        help="Three labels in order: PL, TGglyc, TGacyl")
    p.add_argument('--colors', nargs=3, default=DEFAULT_COLORS,
        help="Three colors for the scatter points")
    p.add_argument('--markers', nargs=3, default=DEFAULT_MARKERS,
        help="Three markers for the scatter points")
    p.add_argument('--title', default=None,
        help="Optional plot title")
    args = p.parse_args()

    d = Path(args.run_dir)
    # adjust these filenames 
    paths = [
        d / 'PLacyl.dat',
        d / 'TGglyc.dat',
        d / 'TGacyl.dat',
    ]
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        p.error(f"Missing files: {', '.join(missing)}")

    plot_defect_data(
        file_paths = paths,
        labels     = args.labels,
        colors     = args.colors,
        markers    = args.markers,
        title      = args.title,
        output_path= args.out
    )


if __name__ == '__main__':
    main()
