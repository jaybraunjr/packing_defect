"""Statistics utilities for packing_defect."""

import os
import argparse
import statistics as stats
from collections import Counter, defaultdict

import numpy as np
from MDAnalysis import Universe

_NEIGHBORS = [
    (dx, dy)
    for dx in (-1, 0, 1)
    for dy in (-1, 0, 1)
    if (dx, dy) != (0, 0)
]


def open_gro(gro_path):
    """
    Open a GRO file and return a set of integer (i, j) grid positions.
    """
    u = Universe(gro_path)
    xy = u.atoms.positions[:, :2]
    ij = np.rint(xy).astype(int)
    return set(map(tuple, ij))


def cluster_positions(positions):
    """
    Return a list of clusters, each as a set of (i, j) tuples.
    """
    visited = set()
    clusters = []
    for pos in positions:
        if pos in visited:
            continue
        stack = [pos]
        visited.add(pos)
        cluster = {pos}
        while stack:
            i0, j0 = stack.pop()
            for di, dj in _NEIGHBORS:
                nb = (i0 + di, j0 + dj)
                if nb in positions and nb not in visited:
                    visited.add(nb)
                    stack.append(nb)
                    cluster.add(nb)
        clusters.append(cluster)
    return clusters


def cluster_sizes(positions):
    """
    Return cluster sizes from a set of positions.
    """
    return [len(c) for c in cluster_positions(positions)]


def count_defect_types(gro_dir, cutoff):
    """
    Count number of clusters ≥ cutoff in each frame GRO file.
    Returns list of (frame_name, count).
    """
    results = []
    for fn in sorted(os.listdir(gro_dir)):
        if not fn.endswith('.gro'):
            continue
        path = os.path.join(gro_dir, fn)
        positions = open_gro(path)
        sizes = cluster_sizes(positions)
        count = sum(1 for s in sizes if s >= cutoff)
        frame_name = fn[:-4]
        results.append((frame_name, count))
    return results


def return_stats(gro_dir, cutoff):
    """
    Return summary statistics (mean, median, stdev) of defect counts.
    """
    results = count_defect_types(gro_dir, cutoff)
    counts = [cnt for (_name, cnt) in results]
    return {
        'mean': stats.mean(counts),
        'median': stats.median(counts),
        'stdev': stats.stdev(counts),
        'per_type': results,
    }


def defect_distribution(gro_dir):
    """
    Returns dict mapping frame_name to Counter of defect sizes.
    """
    distribution = {}
    for fn in sorted(os.listdir(gro_dir)):
        if not fn.endswith('.gro'):
            continue
        path = os.path.join(gro_dir, fn)
        positions = open_gro(path)
        sizes = cluster_sizes(positions)
        distribution[fn[:-4]] = Counter(sizes)
    return distribution


def defect_distro_stats(distros):
    """
    Compute mean and stdev for each defect size across frames.
    """
    size_timeseries = defaultdict(list)
    for distro in distros.values():
        for size, count in distro.items():
            size_timeseries[size].append(count)
    output = {}
    for size, counts in size_timeseries.items():
        output[size] = {
            'mean': stats.mean(counts),
            'stdev': stats.stdev(counts) if len(counts) > 1 else 0.0,
        }
    return output


def average_std_per_frame(distros):
    """
    Return dict mapping frame_name to (mean, stdev) of all defects in that frame.
    """
    result = {}
    for frame, distro in distros.items():
        all_sizes = []
        for size, count in distro.items():
            all_sizes.extend([size] * count)
        if not all_sizes:
            result[frame] = (0.0, 0.0)
        elif len(all_sizes) == 1:
            result[frame] = (all_sizes[0], 0.0)
        else:
            result[frame] = (stats.mean(all_sizes), stats.stdev(all_sizes))
    return result


def global_avg_std_defect_size(distros):
    """
    Return global (mean, stdev) across all frames and sizes.
    """
    all_sizes = []
    for distro in distros.values():
        for size, count in distro.items():
            all_sizes.extend([size] * count)
    if not all_sizes:
        return 0.0, 0.0
    if len(all_sizes) == 1:
        return all_sizes[0], 0.0
    return stats.mean(all_sizes), stats.stdev(all_sizes)


def main():
    parser = argparse.ArgumentParser(
        description="Count defect clusters and compute stats"
    )
    parser.add_argument('gro_dir', help="Directory with .gro files")
    parser.add_argument('-c', '--cutoff', type=int, default=5,
                        help="Minimum cluster size to count")
    parser.add_argument('--summary', action='store_true',
                        help="Print summary stats of counts")
    parser.add_argument('--stats', action='store_true',
                        help="Print per-frame and global defect stats")
    args = parser.parse_args()

    counts = count_defect_types(args.gro_dir, args.cutoff)
    for frame, cnt in counts:
        print(f"{frame:10s} clusters ≥ {args.cutoff:3d}  →  {cnt}")

    if args.summary:
        summary = return_stats(args.gro_dir, args.cutoff)
        print("\nSummary:")
        print(f"  Mean   : {summary['mean']:.2f}")
        print(f"  Median : {summary['median']:.2f}")
        print(f"  Stdev  : {summary['stdev']:.2f}")

    if args.stats:
        distros = defect_distribution(args.gro_dir)
        print("\nAverage defect size per frame:")
        frame_stats = average_std_per_frame(distros)
        for frame, (mean_, std_) in sorted(frame_stats.items()):
            print(f"{frame:20s}  mean: {mean_:.2f}, stdev: {std_:.2f}")

        global_mean, global_std = global_avg_std_defect_size(distros)
        print("\nGlobal average defect size:")
        print(f"  Mean  : {global_mean:.2f}")
        print(f"  Stdev : {global_std:.2f}")

if __name__ == '__main__':
    main()
