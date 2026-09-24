"""Write the cell tables config/cells_{small,paper,full}.tsv.

One row per grid cell. n_chunks splits the trials of a cell into jobs of
similar cost: the cost of a trial grows with n_obs, so larger cells get
more chunks. With obs_per_chunk = 256 and 2000 trials, a chunk of a cell
with n_obs >= 256 takes about 25 min, close to the 30 min per job that
Triton recommends (https://scicomp.aalto.fi/triton/tut/array/). Run from the pipeline folder:
    python workflow/scripts/make_cells.py
"""
import csv
import itertools
import math

# One sub-grid per family. The binomial cells keep beta_t on the logit
# scale, where the same values again run from equal models (0) to clearly
# different ones (1). out_dev shifts the linear predictor of the first
# observation, so its scale is the logit scale too: plogis saturates well
# before 20, and the gaussian levels would be indistinguishable.
GRIDS = {
    # full grid as in simulated/setup_all.py
    "config/cells_full.tsv": dict(
        gaussian=dict(
            n_obs=[16, 32, 64, 128, 256, 512, 1024],
            beta_t=[0.0, 0.05, 0.1, 0.2, 0.5, 1.0],
            out_dev=[0.0, 20.0, 200.0],
            tau2=["none", 1.0],
            obs_per_chunk=256,
        ),
        binomial=dict(
            n_obs=[16, 32, 64, 128, 256, 512, 1024],
            beta_t=[0.0, 0.05, 0.1, 0.2, 0.5, 1.0],
            out_dev=[0.0, 2.0, 5.0],
            tau2=["none"],
            obs_per_chunk=1024,
        ),
        poisson=dict(
            n_obs=[16, 32, 64, 128, 256, 512, 1024],
            beta_t=[0.0, 0.05, 0.1, 0.2, 0.5, 1.0],
            out_dev=[0.0, 1.0, 2.0],
            tau2=["none"],
            obs_per_chunk=256,
        ),
    ),
    # cells of Fig. 5 and Fig. 6 in the paper
    "config/cells_paper.tsv": dict(
        gaussian=dict(
            n_obs=[32, 128, 512],
            beta_t=[0.0, 0.2, 1.0],
            out_dev=[0.0, 20.0],
            tau2=["none"],
            obs_per_chunk=64,
        ),
        binomial=dict(
            n_obs=[32, 128, 512],
            beta_t=[0.0, 0.2, 1.0],
            out_dev=[0.0, 5.0],
            tau2=["none"],
            obs_per_chunk=512,
        ),
        poisson=dict(
            n_obs=[32, 128, 512],
            beta_t=[0.0, 0.2, 1.0],
            out_dev=[0.0, 2.0],
            tau2=["none"],
            obs_per_chunk=128,
        ),
    ),
    # small grid for a laptop test; n_obs = 64 gets 2 chunks
    "config/cells_small.tsv": dict(
        gaussian=dict(
            n_obs=[16, 64],
            beta_t=[0.0, 1.0],
            out_dev=[0.0, 20.0],
            tau2=["none", 1.0],
            obs_per_chunk=32,
        ),
        binomial=dict(
            n_obs=[16, 64],
            beta_t=[0.0, 1.0],
            out_dev=[0.0, 5.0],
            tau2=["none"],
            obs_per_chunk=64,
        ),
        poisson=dict(
            n_obs=[16, 64],
            beta_t=[0.0, 1.0],
            out_dev=[0.0, 2.0],
            tau2=["none"],
            obs_per_chunk=32,
        ),
    ),
}

# Cell names of the gaussian family are unprefixed, so that the cells of
# earlier runs keep their names and their results.
PREFIX = {"gaussian": "", "binomial": "bin_", "poisson": "pois_"}


def fmt(x):
    return x if isinstance(x, str) else f"{x:g}"


for path, families in GRIDS.items():
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            ["cell", "family", "n_obs", "beta_t", "out_dev", "tau2", "n_chunks"]
        )
        for family, grid in families.items():
            for n, b, o, t in itertools.product(
                grid["n_obs"], grid["beta_t"], grid["out_dev"], grid["tau2"]
            ):
                cell = f"{PREFIX[family]}n{n}_b{fmt(b)}_o{fmt(o)}_t{fmt(t)}"
                n_chunks = math.ceil(n / grid["obs_per_chunk"])
                writer.writerow([cell, family, n, fmt(b), fmt(o), fmt(t), n_chunks])
