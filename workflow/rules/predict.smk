import os

# Predict per organs using joined features
rule predict:
    input:
        generic_features = os.path.join(
           BASE, "{run}", "work", "generic_features.parquet"
        ),
        organsp_features = os.path.join(
            BASE, "{run}", "work", "organsp_features",
            "{organ}_features.parquet"
        )
    output:
        file = os.path.join(
            BASE, "{run}", "predictions",
            "TLand_scores.{organ}.tsv.gz"
        )
    log:
        os.path.join(
            BASE, "{run}", "logs",
            "predict.{organ}.log"
        )
    conda:
        "../envs/TLand.yml"
    resources:
        mem_mb=config['memory_predict_organsp_mb']
    shell:
        """
        python workflow/scripts/predict.py \
           --generic_features {input.generic_features} \
           --organsp_features {input.organsp_features} \
           --organ {wildcards.organ} \
           --organ_mapping_json {config[organ_mapping_json]} \
           --organ_list {config[organ_list_path]} \
           --models_path {config[models_path]} \
           --out {output.file} &> {log}
        """