rule summarise:
    input:
        f"{RES}/trials.rds",
    output:
        f"{RES}/moments.rds",
    log:
        f"{RES}/logs/summarise.log",
    conda:
        "../envs/report.yaml"
    params:
        n_bb=2000,
        seed=243682451,
    script:
        "../scripts/summarise.R"


rule plot_calibration:
    input:
        f"{RES}/trials.rds",
    output:
        f"{RES}/figs/calibration_{{family}}_{{measure}}.pdf",
    log:
        f"{RES}/logs/plot_calibration_{{family}}_{{measure}}.log",
    conda:
        "../envs/report.yaml"
    params:
        sel=plot_selection,
    script:
        "../scripts/plot_calibration.R"


rule plot_coverage:
    input:
        f"{RES}/trials.rds",
    output:
        raw=f"{RES}/figs/coverage_absdiff.pdf",
        scaled=f"{RES}/figs/coverage_absdiff_n.pdf",
        cutoffs=f"{RES}/coverage_cutoffs.csv",
        **{
            f"cutoff_{family}": f"{RES}/figs/coverage_cutoff_{family}.pdf"
            for family in FAMILIES
        },
    log:
        f"{RES}/logs/plot_coverage.log",
    conda:
        "../envs/report.yaml"
    params:
        families=FAMILIES,
        n_bins=30,
        target_cov=0.9,
    script:
        "../scripts/plot_coverage.R"


rule plot_joint:
    input:
        f"{RES}/trials.rds",
    output:
        f"{RES}/figs/joint_{{family}}_{{measure}}.pdf",
    log:
        f"{RES}/logs/plot_joint_{{family}}_{{measure}}.log",
    conda:
        "../envs/report.yaml"
    params:
        sel=plot_selection,
    script:
        "../scripts/plot_joint.R"


rule plot_moments:
    input:
        f"{RES}/moments.rds",
    output:
        f"{RES}/figs/moments_{{family}}_{{measure}}.pdf",
    log:
        f"{RES}/logs/plot_moments_{{family}}_{{measure}}.log",
    conda:
        "../envs/report.yaml"
    params:
        sel=plot_selection,
    script:
        "../scripts/plot_moments.R"


rule side_by_side:
    input:
        paper=paper_figure,
        rebuild=f"{RES}/figs/{{fig}}_gaussian_elpd.pdf",
    output:
        f"{RES}/figs/side_by_side_{{fig}}.pdf",
    log:
        f"{RES}/logs/side_by_side_{{fig}}.log",
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/side_by_side.R"


rule plot_pointwise:
    input:
        a=cell_files("{res}/{cell}/score/A_{chunk}.rds"),
        b=cell_files("{res}/{cell}/score/B_{chunk}.rds"),
        compare=cell_files("{res}/{cell}/compare/{chunk}.rds"),
    output:
        f"{RES}/figs/pointwise_{{cell}}.pdf",
    log:
        f"{RES}/logs/plot_pointwise_{{cell}}.log",
    conda:
        "../envs/report.yaml"
    params:
        cell=cell_params,
        measures=config["measures"],
    script:
        "../scripts/plot_pointwise.R"


rule plot_var_ratio:
    input:
        f"{RES}/trials.rds",
    output:
        f"{RES}/figs/var_ratio_{{family}}_{{measure}}.pdf",
    log:
        f"{RES}/logs/plot_var_ratio_{{family}}_{{measure}}.log",
    conda:
        "../envs/report.yaml"
    params:
        sel=plot_selection,
    script:
        "../scripts/plot_var_ratio.R"


rule plot_err:
    input:
        f"{RES}/trials.rds",
    output:
        f"{RES}/figs/err_{{family}}_{{measure}}.pdf",
    log:
        f"{RES}/logs/plot_err_{{family}}_{{measure}}.log",
    conda:
        "../envs/report.yaml"
    params:
        sel=plot_selection,
        signed=False,
    script:
        "../scripts/plot_err.R"


rule plot_errdirection:
    input:
        f"{RES}/trials.rds",
    output:
        f"{RES}/figs/errdirection_{{family}}_{{measure}}.pdf",
    log:
        f"{RES}/logs/plot_errdirection_{{family}}_{{measure}}.log",
    conda:
        "../envs/report.yaml"
    params:
        sel=plot_selection,
        signed=True,
    script:
        "../scripts/plot_err.R"


rule plot_overview:
    input:
        trials=f"{RES}/trials.rds",
        moments=f"{RES}/moments.rds",
    output:
        f"{RES}/figs/overview_{{family}}_{{measure}}.pdf",
    log:
        f"{RES}/logs/plot_overview_{{family}}_{{measure}}.log",
    conda:
        "../envs/report.yaml"
    params:
        sel=plot_selection,
    script:
        "../scripts/plot_overview.R"
