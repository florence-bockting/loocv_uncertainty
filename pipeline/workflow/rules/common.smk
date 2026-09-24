import csv


def read_cells(path):
    """Read the cell table: one row per grid cell."""
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return {
        row["cell"]: dict(
            family=row.get("family", "gaussian"),
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
FAMILIES = sorted({cell["family"] for cell in CELLS.values()})
# Mirrors MEASURES_BINARY in workflow/scripts/lib/measures.R.
MEASURES_BINARY = ["brier", "acc", "bacc"]


def family_measures(family):
    """The configured measures that a family supports."""
    return [
        m
        for m in config["measures"]
        if family == "binomial" or m not in MEASURES_BINARY
    ]


wildcard_constraints:
    cell="|".join(CELLS),
    model="|".join(MODELS),
    family="|".join(FAMILIES),
    chunk=r"\d+",
    measure=r"[a-z0-9]+",
    fig=r"[a-z_]+",


def cell_params(wildcards):
    return CELLS[wildcards.cell]


def chunk_trials(wildcards):
    """First and last trial (1-based, inclusive) of a chunk."""
    n_trial = config["n_trial"]
    n_chunks = CELLS[wildcards.cell]["n_chunks"]
    k = int(wildcards.chunk)
    return [(k - 1) * n_trial // n_chunks + 1, k * n_trial // n_chunks]


# "overview" is the one-page summary; the others have a paper counterpart.
FIGS = [
    "overview",
    "calibration",
    "joint",
    "moments",
    "var_ratio",
    "err",
    "errdirection",
]


def plot_selection(wildcards):
    """The cells a figure shows. `plot_<family>` in the config overrides
    single entries of `plot` for that family."""
    sel = dict(config["plot"])
    sel.update(config.get(f"plot_{wildcards.family}", {}))
    return sel


def plot_cells():
    """The cells of the plot selection, for the figures of a single cell."""
    cells = []
    for cell, params in CELLS.items():
        sel = dict(config["plot"])
        sel.update(config.get(f"plot_{params['family']}", {}))
        if (
            params["n_obs"] in sel["n_obs"]
            and params["beta_t"] in [float(b) for b in sel["beta_t"]]
            and params["out_dev"] in [float(o) for o in sel["out_dev"]]
            and params["tau2"] == str(sel["tau2"])
        ):
            cells.append(cell)
    return cells


def cell_files(pattern):
    """`pattern` formatted with the cell and every chunk of that cell."""

    def files(wildcards):
        n_chunks = CELLS[wildcards.cell]["n_chunks"]
        return [
            pattern.format(res=RES, cell=wildcards.cell, chunk=chunk)
            for chunk in range(1, n_chunks + 1)
        ]

    return files


def paper_figure(wildcards):
    return config["paper_figs"][wildcards.fig]


def report_files():
    """Figures for every measure, plus the side-by-side pages if the config
    names the paper figures."""
    files = [
        f"{RES}/figs/{fig}_{family}_{measure}.pdf"
        for fig in FIGS
        for family in FAMILIES
        for measure in family_measures(family)
    ]
    files += [f"{RES}/figs/pointwise_{cell}.pdf" for cell in plot_cells()]
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
