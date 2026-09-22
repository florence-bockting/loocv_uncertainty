import csv


def read_cells(path):
    """Read the cell table: one row per grid cell."""
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return {
        row["cell"]: dict(
            n_obs=int(row["n_obs"]),
            beta_t=float(row["beta_t"]),
            out_dev=float(row["out_dev"]),
            tau2=row["tau2"],
            n_chunks=int(row["n_chunks"]),
        )
        for row in rows
    }


CELLS = read_cells(config["cells"])
RES = config["results"]
N_OBS_MAX = max(cell["n_obs"] for cell in CELLS.values())
MODELS = ["A", "B"]


wildcard_constraints:
    cell="|".join(CELLS),
    model="|".join(MODELS),
    chunk=r"\d+",
    measure=r"[a-z0-9]+",
    fig=r"[a-z]+",


def cell_params(wildcards):
    return CELLS[wildcards.cell]


def chunk_trials(wildcards):
    """First and last trial (1-based, inclusive) of a chunk."""
    n_trial = config["n_trial"]
    n_chunks = CELLS[wildcards.cell]["n_chunks"]
    k = int(wildcards.chunk)
    return [(k - 1) * n_trial // n_chunks + 1, k * n_trial // n_chunks]


FIGS = ["calibration", "joint", "moments"]


def paper_figure(wildcards):
    return config["paper_figs"][wildcards.fig]


def report_files():
    """Figures for every measure, plus the side-by-side pages if the config
    names the paper figures."""
    files = [
        f"{RES}/figs/{fig}_{measure}.pdf"
        for fig in FIGS
        for measure in config["measures"]
    ]
    files += [
        f"{RES}/figs/side_by_side_{fig}.pdf" for fig in config.get("paper_figs", {})
    ]
    return files


def all_compare_files(wildcards):
    return [
        f"{RES}/{cell}/compare/{chunk}.rds"
        for cell, params in CELLS.items()
        for chunk in range(1, params["n_chunks"] + 1)
    ]
