

import numpy as np
import os
import argparse
from MDAnalysis import Universe
import statistics
from collections import Counter
from collections import defaultdict


_NEIGHBORS = [(dx,dy)
              for dx in (-1,0,1)
              for dy in (-1,0,1)
            if not (dx==0 and dy==0)]

def open_gro(gro_path):

    u = Universe(gro_path)
    xy = u.atoms.positions[:, :2]
    ij = np.rint(xy).astype(int)
    return set(map(tuple, ij))



def cluster_positions(positions):
    """Return a list of clusters, each as a set of (i,j) grid tuples."""
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
    """Legacy: only return lengths."""
    return [len(c) for c in cluster_positions(positions)]



def count_defect_types(gro_dir, cutoff):

    results = []
    for fn in sorted(os.listdir(gro_dir)):
        if not fn.endswith('.gro'):
            continue
        path = os.path.join(gro_dir, fn)
        positions = open_gro(path)
        sizes = cluster_sizes(positions)
        # count how many clusters meet the cutoff
        num = sum(1 for s in sizes if s >= cutoff)
        results.append((fn[:-4], num))
    return results




def return_stats(gro_dir, cutoff):

    results = count_defect_types(gro_dir, cutoff)
    counts = [cnt for (_name, cnt) in results]

    mean = statistics.mean(counts)
    median = statistics.median(counts)
    stdev = statistics.stdev(counts)

    # print("mean:",mean, "median:",median, "stdev:",stdev )

    return {"mean": mean,
            "median": median,
            "stdev": stdev,
            "per_type": results
            }



def defect_distribution(gro_dir):
    """ 
    Returns a dict where each key is a frame name and each value
    is a counter of defect sizes observed in that frame
    """
    frame_distro = {}

    for fn in sorted(os.listdir(gro_dir)):
        if not fn.endswith('.gro'):
            continue
        path = os.path.join(gro_dir, fn)
        positions = open_gro(path)
        sizes = cluster_sizes(positions)
        size_count = Counter(sizes)
        key = fn[:-4]
        frame_distro[key] = size_count
    return(frame_distro)


def defect_distro_stats(defect_distribution):

    size_timeseries = defaultdict(list)
    for frame_name, size_counter in defect_distribution.items():
        for size, count in size_counter.items():
            size_timeseries[size].append(count)

    output = {}
    for size, counts in size_timeseries.items():
        output[size] = {
            'mean': statistics.mean(counts),
            'stdev': statistics.stdev(counts) if len(counts) > 1 else 0.0
        }

    return output


def average_std_per_frame(defect_distribution):
    """
    Returns a dict: {frame_name: (mean, stdev)} for defect sizes in that frame
    """
    from statistics import mean, stdev

    result = {}
    for frame_name, size_counter in defect_distribution.items():
        all_sizes = []
        for size, count in size_counter.items():
            all_sizes.extend([size] * count)

        if len(all_sizes) == 0:
            result[frame_name] = (0.0, 0.0)
        elif len(all_sizes) == 1:
            result[frame_name] = (all_sizes[0], 0.0)
        else:
            result[frame_name] = (mean(all_sizes), stdev(all_sizes))
    return result




def global_avg_std_defect_size(defect_distribution):
    """
    Returns global (mean, stdev) across all frames and all sizes
    """
    from statistics import mean, stdev

    all_sizes = []
    for size_counter in defect_distribution.values():
        for size, count in size_counter.items():
            all_sizes.extend([size] * count)

    if len(all_sizes) == 0:
        return (0.0, 0.0)
    elif len(all_sizes) == 1:
        return (all_sizes[0], 0.0)
    else:
        return (mean(all_sizes), stdev(all_sizes))





def main():
    p = argparse.ArgumentParser(
        description="Count how many defect clusters exceed a size cutoff"
    )
    p.add_argument('gro_dir',
                   help="Directory containing .gro files from run_defect")
    p.add_argument('-c','--cutoff', type=int, default=5,
                   help="Minimum cluster size (in grid cells) to count")
    p.add_argument('--summary', action='store_true',
                   help="Print mean, median, stdev of defect counts")
    p.add_argument('--stats', action='store_true',
                   help="Calculate average and std of defect sizes per frame and globally")
    args = p.parse_args()

    stats = count_defect_types(args.gro_dir, args.cutoff)
    for name, cnt in stats:
        print(f"{name:10s} clusters ≥ {args.cutoff:3d}  →  {cnt}")

    if args.summary:
        summary = return_stats(args.gro_dir, args.cutoff)
        print("\nSummary:")
        print(f"  Mean   : {summary['mean']:.2f}")
        print(f"  Median : {summary['median']:.2f}")
        print(f"  Stdev  : {summary['stdev']:.2f}")

    if args.stats:
        if args.stats:
            defect_distro = defect_distribution(args.gro_dir)

            print("\nAverage defect size per frame:")
            frame_avgs = average_std_per_frame(defect_distro)
            for frame, (avg, std) in sorted(frame_avgs.items()):
                print(f"{frame:20s}  mean: {avg:.2f}, stdev: {std:.2f}")

        global_avg, global_std = global_avg_std_defect_size(defect_distro)
        print("\nGlobal average defect size:")
        print(f"  Mean  : {global_avg:.2f}")
        print(f"  Stdev : {global_std:.2f}")


if __name__ == '__main__':
    main()
