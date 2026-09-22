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
        f"{RES}/figs/calibration_{{measure}}.pdf",
    log:
        f"{RES}/logs/plot_calibration_{{measure}}.log",
    conda:
        "../envs/report.yaml"
    params:
        sel=config["plot"],
    script:
        "../scripts/plot_calibration.R"


rule plot_joint:
    input:
        f"{RES}/trials.rds",
    output:
        f"{RES}/figs/joint_{{measure}}.pdf",
    log:
        f"{RES}/logs/plot_joint_{{measure}}.log",
    conda:
        "../envs/report.yaml"
    params:
        sel=config["plot"],
    script:
        "../scripts/plot_joint.R"


rule plot_moments:
    input:
        f"{RES}/moments.rds",
    output:
        f"{RES}/figs/moments_{{measure}}.pdf",
    log:
        f"{RES}/logs/plot_moments_{{measure}}.log",
    conda:
        "../envs/report.yaml"
    params:
        sel=config["plot"],
    script:
        "../scripts/plot_moments.R"


rule side_by_side:
    input:
        paper=paper_figure,
        rebuild=f"{RES}/figs/{{fig}}_elpd.pdf",
    output:
        f"{RES}/figs/side_by_side_{{fig}}.pdf",
    log:
        f"{RES}/logs/side_by_side_{{fig}}.log",
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/side_by_side.R"
