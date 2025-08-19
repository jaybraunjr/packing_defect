import os
import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np

class LeafletVisualizer:

    def __init__(self, base_dir, lipid, leaflet="up", method="box", dpi=150):
        self.base_dir = base_dir
        self.lipid = lipid
        self.leaflet = leaflet
        self.method = method
        self.dpi = dpi

    def _z_mid(self, u):
        if self.method =='box':
            return float(u.dimensions[2]) * 0.5
        elif self.method =='data':
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


    def _leaflet_xy(self, u):
        mask = self._leaflet_mask(u)
        pos = u.atoms.positions[mask]
        return pos[:, 0], pos[:, 1]


    def render_frame(self, frame):
        gro = os.path.join(self.base_dir, self.lipid, f"{self.lipid}_frame_{frame}.gro")
        if not os.path.exists(gro):
            return None

        u = mda.Universe(gro)
        fig, ax = plt.subplots(figsize=(6, 6))
        x, y = self._leaflet_xy(u)
        ax.scatter(x, y, s=3, alpha=0.8)
        Lx, Ly = u.dimensions[:2]
        if Lx and Ly:
            ax.set_xlim(0, Lx)
            ax.set_ylim(0, Ly)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x (Å)")
        ax.set_ylabel("y (Å)")
        ax.set_title(f"{self.lipid} — frame {frame} — {self.leaflet}")

        outdir = os.path.join(self.base_dir, f"frames_{self.lipid}_{self.leaflet}")
        os.makedirs(outdir, exist_ok=True)
        outpng = os.path.join(outdir, f"{self.lipid}_{self.leaflet}_frame_{frame:04d}.png")

        plt.savefig(outpng, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)
        return outpng



    def render_frames(self, start, stop):
        outputs = []
        for f in range(start, stop):
            out = self.render_frame(f)
            if out:
                outputs.append(out)
        return outputs