import os
import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np


class LeafletVisualizer:
    """
    Visualize lipid defect positions separated by leaflet, with optional protein backbone overlay.
    """

    def __init__(self, base_dir, lipid, leaflet="up", method="box", dpi=150,
                 defect_color="blue", protein_color="red"):
        self.base_dir = base_dir
        self.lipid = lipid
        self.leaflet = leaflet
        self.method = method
        self.dpi = dpi
        self.defect_color = defect_color
        self.protein_color = protein_color

    # ---------------- Internal helpers ----------------
    def _z_mid(self, u):
        if self.method == "box":
            return float(u.dimensions[2]) * 0.5
        elif self.method == "data":
            z = u.atoms.positions[:, 2]
            return 0.5 * (float(z.min()) + float(z.max()))
        else:
            raise ValueError("method must be 'box' or 'data'")

    def _leaflet_mask(self, u):
        mid = self._z_mid(u)
        z = u.atoms.positions[:, 2]
        if self.leaflet == "up":
            return z > mid
        elif self.leaflet == "dw":
            return z <= mid
        else:
            raise ValueError("leaflet must be 'up' or 'dw'")



    def _defect_xy(self, u):
        """Defect atoms = all non-protein atoms split by leaflet."""
        sel = u.select_atoms("not protein")
        mid = self._z_mid(u)
        z = sel.positions[:, 2]
        mask = (z > mid) if self.leaflet == "up" else (z <= mid)
        pos = sel.positions[mask]
        return pos[:, 0], pos[:, 1]

    def _protein_xy(self, u):
        """Protein backbone atoms restricted to the chosen leaflet."""
        sel = u.select_atoms("protein and backbone")
        if sel.n_atoms == 0:
            return np.array([]), np.array([])

        mid = self._z_mid(u)
        z = sel.positions[:, 2]
        mask = (z > mid) if self.leaflet == "up" else (z <= mid)
        pos = sel.positions[mask]
        return pos[:, 0], pos[:, 1]


    # ---------------- Plotting ----------------
    def _plot(self, u, frame, show=True, savepath=None):
        fig, ax = plt.subplots(figsize=(6, 6))

        # defects
        x, y = self._defect_xy(u)
        ax.scatter(x, y, s=3, alpha=0.8, c=self.defect_color, label="Defects")

        # protein backbone
        x, y = self._protein_xy(u)
        if x.size > 0:
            ax.plot(x, y, "o-", color=self.protein_color,
                    markersize=2, linewidth=0.5, label="Protein backbone")

        Lx, Ly = u.dimensions[:2]
        if Lx and Ly:
            ax.set_xlim(0, Lx)
            ax.set_ylim(0, Ly)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x (Å)")
        ax.set_ylabel("y (Å)")
        ax.set_title(f"{self.lipid} — frame {frame} — {self.leaflet}")
        ax.legend()

        if savepath:
            plt.savefig(savepath, dpi=self.dpi, bbox_inches="tight")
            plt.close(fig)
            return savepath

        if show:
            plt.show()

        return fig, ax

    def render_frame(self, frame):
        """Save a single frame to PNG."""
        gro = os.path.join(self.base_dir, self.lipid, f"{self.lipid}_frame_{frame}.gro")
        if not os.path.exists(gro):
            return None

        u = mda.Universe(gro)
        outdir = os.path.join(self.base_dir, f"frames_{self.lipid}_{self.leaflet}")
        os.makedirs(outdir, exist_ok=True)
        outpng = os.path.join(outdir, f"{self.lipid}_{self.leaflet}_frame_{frame:04d}.png")
        return self._plot(u, frame, show=False, savepath=outpng)

    def render_frames(self, start, stop):
        """Save a range of frames to PNGs."""
        outputs = []
        for f in range(start, stop):
            out = self.render_frame(f)
            if out:
                outputs.append(out)
        return outputs

    def plot_frame(self, frame, show=True):
        """Plot a single frame inline in a notebook."""
        gro = os.path.join(self.base_dir, self.lipid, f"{self.lipid}_frame_{frame}.gro")
        if not os.path.exists(gro):
            return None, None

        u = mda.Universe(gro)
        return self._plot(u, frame, show=show)
