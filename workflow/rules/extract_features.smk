import os
from snakemake.io import directory

# 1) Prep input bed file for RegulomeDB query from vcf-like file (chr, pos, id, ref, alt)
rule prep_input_bed:
    input:
        vcf=lambda wc: RUN_PARAMS[wc.run]["VCF"]
    output:
        bed=os.path.join(
            BASE, "{run}", "work", "reg_query_input.bed"
        )
    shell:
        """
        awk -F"\t" 'BEGIN{{OFS="\t"}} {{print $1, $2-1, $2}}' {input.vcf} | sort | uniq > {output.bed}
        """

# 2a) Query variants from the database (requires running from specific directory)
rule query_variants:
    input:
        bed=os.path.join(
            BASE, "{run}", "work", "reg_query_input.bed"
        )
    output:
        jsonl=os.path.join(
            BASE, "{run}", "work", "regdb_query_output.jsonl"
        )
    log:
        os.path.join(
            BASE, "{run}", "logs", "reg_query.log"
        )
    conda:
        "../envs/gds.yml"
    resources:
        queries=1
    shell:
        """
        cd {config[gds_dir]}
        python -m utils.regulome_search_TLand \
            -f {input.bed} \
            --assembly GRCh38 \
            --peaks 1> {output.jsonl} 2> {log}
        """

# 2b) Run Sei on input VCF
rule run_sei_vep:
    input:
        vcf=lambda wc: RUN_PARAMS[wc.run]["VCF"]
    output:
        outdir=temp(directory(os.path.join(
            BASE, "{run}", "work", "sei", "sei_output"
        )))
    log:
        os.path.join(
            BASE, "{run}", "logs", "run_sei_pipeline.log"
        )
    conda:
        "../envs/sei.yml"
    resources:
        gpu=config["gpu"]
    shell:
        """
        cd {config[sei_dir]}
        sh run_pipeline.sh {input.vcf} hg38 {output.outdir} --cuda &> {log}
        """

rule run_sei_seq_class:
    input:
        vcf=lambda wc: RUN_PARAMS[wc.run]["VCF"],
        vep_outdir=os.path.join(
            BASE, "{run}", "work", "sei", "sei_output"
        ),
        sei_model_dir=config["sei_dir"]+"/model"
    output:
        features=os.path.join(
            BASE, "{run}", "work", "sei", "sei_final_output", "sei_features.tsv"
        )
    log:
        os.path.join(
            BASE, "{run}", "logs", "run_seq_class.log"
        )
    conda:
        "../envs/sei.yml"
    shell:
        """
        python workflow/scripts/run_seq_class.py -s {input.vep_outdir} -i {input.vcf} -m {input.sei_model_dir} -o {output.features} &> {log}
        """

# 2) Extract generic features
rule extract_generic:
    input:
        vcf=lambda wc: RUN_PARAMS[wc.run]["VCF"],
        jsonl=os.path.join(
            BASE, "{run}", "work", "regdb_query_output.jsonl"
        ),
        sei_features=os.path.join(
            BASE, "{run}", "work", "sei", "sei_final_output", "sei_features.tsv"
        )
    output:
        parquet=os.path.join(
            BASE, "{run}", "work", "generic_features.parquet"
        )
    log:
        os.path.join(
            BASE, "{run}", "logs", "extract_generic_features.log"
        )
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
        jsonl=os.path.join(
            BASE, "{run}", "work", "regdb_query_output.jsonl"
        ),
        total_num_path=config["total_num_path"]
    output:
        parquet=os.path.join(
            BASE, "{run}", "work", "organsp_features", "{organ}_features.parquet"
        )
    log:
        os.path.join(
            BASE, "{run}", "logs", "extract_{organ}_features.log"
        )
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
           --organ_mapping_json {config[organ_mapping_json]} \
           --out {output.parquet} &> {log}
        """