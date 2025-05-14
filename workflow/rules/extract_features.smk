import os
from snakemake.io import directory

configfile: "config/config.yaml"

# Available organs from config
ORGANS = config["organs"]

# Directories
WORK = os.path.abspath(config["work_dir"]) # Intermediate files
OUT = config["out_dir"] # Final prediction files
LOGS = os.path.abspath(config["log_dir"])

# 1) Prep input bed file for RegulomeDB query from vcf-like file (chr, pos, id, ref, alt)
rule prep_input_bed:
    input:
        vcf=config["input_vcf"]
    output:
        bed=WORK+"/reg_query_input.bed"
    shell:
        """
        awk -F"\t" 'BEGIN{{OFS="\t"}} {{print $1, $2-1, $2}}' {input.vcf} > {output.bed}
        """

# 2a) Query variants from the database (requires running from specific directory)
rule query_variants:
    input:
        bed=WORK+"/reg_query_input.bed"
    output:
        jsonl=WORK+"/regdb_query_output.jsonl"
    log:
        LOGS+"/reg_query.log"
    conda:
        "../envs/gds.yml"
    shell:
        """
        cd {config[gds_dir]}
        python -m utils.regulome_search_organ \
            -f {input.bed} \
            --assembly GRCh38 \
            --peaks 1> {output.jsonl} 2> {log}
        """

# 2b) Run Sei on input VCF
rule run_sei_vep:
    input:
        vcf=lambda wildcards: os.path.abspath(config["input_vcf"])
    output:
        outdir=directory(WORK+"/sei/sei_output")
    log:
        LOGS+"/run_sei_pipeline.log"
    conda:
        "../envs/sei.yml"
    shell:
        """
        cd {config[sei_dir]}
        sh run_pipeline.sh {input.vcf} hg38 {output.outdir} --cuda &> {log}
        """

rule run_sei_seq_class:
    input:
        vcf=config["input_vcf"],
        vep_outdir=WORK+"/sei/sei_output",
        sei_model_dir=config["sei_dir"]+"/model"
    output:
        features=WORK+"/sei/sei_final_output/sei_features.tsv"
    log:
        LOGS+"/run_seq_class.log"
    conda:
        "../envs/sei.yml"
    shell:
        """
        python workflow/scripts/run_seq_class.py -s {input.vep_outdir} -i {input.vcf} -m {input.sei_model_dir} -o {output.features} &> {log}
        """

# 2) Extract generic features
rule extract_generic:
    input:
        vcf=config["input_vcf"],
        jsonl=WORK+"/regdb_query_output.jsonl",
        sei_features=WORK+"/sei/sei_final_output/sei_features.tsv"
    output:
        parquet=WORK+"/generic_features.parquet"
    log:
        LOGS+"/extract_generic_features.log"
    conda:
        "../envs/TLand.yml"
    resources:
        mem_mb=config["memory_extract_generic_mb"]
    shell:
        """
        python workflow/scripts/extract_generic_features.py \
            --input_vcf {input.vcf} \
            --input_jsonl {input.jsonl} \
            --input_sei {input.sei_features} \
            --dnase_sig_path {config[dnase_sig_path]} \
            --chip_sig_path {config[chip_sig_path]} \
            --out {output.parquet} &> {log}
        """
# 3) Extract organ-specific features in parallel
rule extract_organsp_features:
    input:
        jsonl=WORK+"/regdb_query_output.jsonl",
        total_num_path=config["total_num_path"]
    output:
        parquet=WORK+"/organsp_features/{organ}_features.parquet"
    log:
        LOGS+"/extract_{organ}_features.log"
    conda:
        "../envs/TLand.yml"
    resources:
        mem_mb=config["memory_extract_organsp_mb"]
    shell:
        """
        python workflow/scripts/extract_organsp_features.py \
           --input_jsonl {input.jsonl} \
           --organ {wildcards.organ} \
           --organsp_dnase_sig_path {config[organsp_dnase_sig_path]} \
           --total_num_path {input.total_num_path} \
           --out {output.parquet} &> {log}
        """