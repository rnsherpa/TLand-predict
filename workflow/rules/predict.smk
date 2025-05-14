import os

configfile: "config/config.yaml"

# Available organs from config
ORGANS = config["organs"]

# Directories
WORK = os.path.abspath(config["work_dir"]) # Intermediate files
OUT = config["out_dir"] # Final prediction files
LOGS = os.path.abspath(config["log_dir"])

# Predict per organs using joined features
rule predict:
    input:
        generic_features=WORK + "/generic_features.parquet",
        organsp_features=WORK + "/organsp_features/{organ}_features.parquet" 
    output:
        file=OUT + "/TLand_scores.{organ}.tsv.gz"
    log:
        LOGS + "/predict.{organ}.log"
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
           --organ_list {config[organ_list_path]} \
           --models_path {config[models_path]} \
           --out {output.file} &> {log}
        """