rule compare_models:
    input:
        a=grid_files("{res}/{cell}/score/A_{chunk}.rds"),
        b=grid_files("{res}/{cell}/score/B_{chunk}.rds"),
        train=[f"{RES}/{cell}/train.rds" for cell in CELLS],
        info=[f"{RES}/{cell}/data_info.rds" for cell in CELLS],
    output:
        f"{RES}/trials.rds",
    log:
        f"{RES}/logs/compare_models.log",
    conda:
        "../envs/r.yaml"
    params:
        cells=CELLS,
        measures=config["measures"],
    script:
        "../scripts/compare_models.R"
