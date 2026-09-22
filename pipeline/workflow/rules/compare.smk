rule compare:
    input:
        a=f"{RES}/{{cell}}/score/A_{{chunk}}.rds",
        b=f"{RES}/{{cell}}/score/B_{{chunk}}.rds",
        train=f"{RES}/{{cell}}/train.rds",
        info=f"{RES}/{{cell}}/data_info.rds",
    output:
        f"{RES}/{{cell}}/compare/{{chunk}}.rds",
    log:
        f"{RES}/logs/{{cell}}/compare_{{chunk}}.log",
    group:
        "chunk"
    conda:
        "../envs/r.yaml"
    params:
        cell=cell_params,
        measures=config["measures"],
    script:
        "../scripts/compare.R"


rule combine:
    input:
        all_compare_files,
    output:
        f"{RES}/trials.rds",
    log:
        f"{RES}/logs/combine.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/combine.R"
