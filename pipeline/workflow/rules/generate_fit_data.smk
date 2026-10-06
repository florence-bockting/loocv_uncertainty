rule generate_data:
    output:
        train=f"{RES}/{{cell}}/train.rds",
        test=temp(f"{RES}/{{cell}}/test.rds"),
        info=f"{RES}/{{cell}}/data_info.rds",
    log:
        f"{RES}/logs/{{cell}}/data.log",
    conda:
        "../envs/r.yaml"
    params:
        cell=cell_params,
        n_sets_train=config["n_sets_train"],
        n_sets_test_gaussian=config["n_sets_test_gaussian"],
        n_sets_test_binomial=config.get("n_sets_test_binomial", 250),
        n_sets_test_poisson=config.get("n_sets_test_poisson", 50),
        n_obs_max=N_OBS_MAX,
        seed=config["seed"],
    script:
        "../scripts/generate_data.R"


rule fit_posterior_and_loo:
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
