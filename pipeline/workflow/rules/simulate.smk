rule data:
    output:
        train=temp(f"{RES}/{{cell}}/train.rds"),
        test=temp(f"{RES}/{{cell}}/test.rds"),
        info=f"{RES}/{{cell}}/data_info.rds",
    log:
        f"{RES}/logs/{{cell}}/data.log",
    conda:
        "../envs/r.yaml"
    params:
        cell=cell_params,
        n_trial=config["n_trial"],
        n_test=config["n_test"],
        n_obs_max=N_OBS_MAX,
        seed=config["seed"],
    script:
        "../scripts/data.R"


rule fit:
    input:
        train=f"{RES}/{{cell}}/train.rds",
    output:
        f"{RES}/{{cell}}/fit/{{model}}_{{chunk}}.rds",
    log:
        f"{RES}/logs/{{cell}}/fit_{{model}}_{{chunk}}.log",
    group:
        "chunk"
    conda:
        "../envs/r.yaml"
    params:
        cell=cell_params,
        trials=chunk_trials,
    script:
        "../scripts/fit.R"
