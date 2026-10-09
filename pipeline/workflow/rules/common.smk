import csv
import math


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
            # Inf: normal covariate; finite: t covariate (gaussian only)
            x_df=float(row.get("x_df", "Inf")),
            n_chunks=int(row["n_chunks"]),
        )
        for row in rows
    }


CELLS = read_cells(config["cells"])
RES = config["results"]
N_OBS_MAX = max(cell["n_obs"] for cell in CELLS.values())
MODELS = ["A", "B"]
FAMILIES = sorted({cell["family"] for cell in CELLS.values()})
# Degrees of freedom of the t covariate in the grid; the normal is not listed.
X_DFS = sorted(
    {cell["x_df"] for cell in CELLS.values() if math.isfinite(cell["x_df"])}
)
# Mirror MEASURES_BINARY and MEASURES_GAUSSIAN in
# workflow/scripts/lib/measures.R.
MEASURES_BINARY = ["brier", "acc", "bacc"]
MEASURES_GAUSSIAN = ["rps", "srps"]


def family_measures(family):
    """The configured measures that a family supports."""
    return [
        m
        for m in config["measures"]
        if (family == "binomial" or m not in MEASURES_BINARY)
        and (family == "gaussian" or m not in MEASURES_GAUSSIAN)
    ]


wildcard_constraints:
    cell="|".join(CELLS),
    model="|".join(MODELS),
    family="|".join(FAMILIES),
    chunk=r"\d+",
    measure=r"[a-z0-9]+",
    fig=r"[a-z_]+",
    # empty for the normal covariate, "_t<df>" for a t covariate
    xdf=r"(_t[0-9.]+)?",


def cell_params(wildcards):
    return CELLS[wildcards.cell]


def chunk_trials(wildcards):
    """First and last trial (1-based, inclusive) of a chunk."""
    n_sets_train = config["n_sets_train"]
    n_chunks = CELLS[wildcards.cell]["n_chunks"]
    k = int(wildcards.chunk)
    return [(k - 1) * n_sets_train // n_chunks + 1, k * n_sets_train // n_chunks]


# "overview" is the one-page summary; the others have a paper counterpart.
PAPER_FIGS = [
    "calibration",
    "joint",
    "moments",
    "var_ratio",
    "err",
    "errdirection",
]
FIGS = ["overview"] + PAPER_FIGS


def x_df_param(wildcards):
    """x_df of a figure: Inf for the normal covariate."""
    return float(wildcards.xdf[2:]) if wildcards.xdf else math.inf


def xdf_suffix(x_df):
    return f"_t{x_df:g}"


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
        if math.isfinite(params["x_df"]):
            continue
        sel = dict(config["plot"])
        sel.update(config.get(f"plot_{params['family']}", {}))
        if (
            params["n_obs"] in sel["n_obs"]
            and params["beta_t"] in [float(b) for b in sel["beta_t"]]
            and params["out_dev"] in [float(o) for o in sel["out_dev"]]
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
    names the paper figures. A t covariate gets its own paper figures; all
    other figures show the normal covariate only."""
    files = [
        f"{RES}/figs/{fig}_{family}_{measure}.pdf"
        for fig in FIGS
        for family in FAMILIES
        for measure in family_measures(family)
    ]
    files += [
        f"{RES}/figs/{fig}_gaussian_{measure}{xdf_suffix(x_df)}.pdf"
        for fig in PAPER_FIGS
        for measure in family_measures("gaussian")
        for x_df in X_DFS
    ]
    files += [f"{RES}/figs/pointwise_{cell}.pdf" for cell in plot_cells()]
    files += [
        f"{RES}/figs/coverage_absdiff.pdf",
        f"{RES}/figs/coverage_absdiff_n.pdf",
        f"{RES}/coverage_cutoffs.csv",
    ]
    files += [f"{RES}/figs/coverage_cutoff_{family}.pdf" for family in FAMILIES]
    files += [
        f"{RES}/figs/side_by_side_{fig}.pdf" for fig in config.get("paper_figs", {})
    ]
    return files


def grid_files(pattern):
    """`pattern` formatted with every cell and chunk, in the order of CELLS."""
    return [
        pattern.format(res=RES, cell=cell, chunk=chunk)
        for cell, params in CELLS.items()
        for chunk in range(1, params["n_chunks"] + 1)
    ]
