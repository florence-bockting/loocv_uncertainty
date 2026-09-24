rule score:
    input:
        fit=f"{RES}/{{cell}}/fit/{{model}}_{{chunk}}.rds",
        train=f"{RES}/{{cell}}/train.rds",
        test=f"{RES}/{{cell}}/test.rds",
    output:
        f"{RES}/{{cell}}/score/{{model}}_{{chunk}}.rds",
    log:
        f"{RES}/logs/{{cell}}/score_{{model}}_{{chunk}}.log",
    group:
        "chunk"
    conda:
        "../envs/r.yaml"
    params:
        cell=cell_params,
        scores=config["scores"],
    script:
        "../scripts/score.R"
