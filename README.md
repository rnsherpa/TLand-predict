# TLand-predict

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)

A Snakemake pipeline to create input feature tables and run predictions with TLand

## Usage

### Input data

Input variant lists should be in a VCF-like tab-delimited file with no headers and columns representing:

1. Chromosome in 'chrN' format
2. 1-based position of variant
3. Variant identifier (e.g. rsID, SPDI, etc.)
4. Reference allele
5. Alternate allele

An example file can be found at `example/example.vcf`

### Runs table

The runs table is a tab-delimited file with three columns: `run_id`, `organs`, and `input_vcf`. An example file can be found at `example/example_runs_table.tsv`

The contents of the run table are described below:

| column name        | details |
| ------------------ | ------- |
| run_id             | The identifier of the run |
| organs             | Desired organs for TLand predictions. Can be any combination of the 51 organs listed in `resources/organ_list.txt` (make sure to use underscores instead of spaces for multi-word organs). Multiple organs must be separated using a semicolon ';'. All 51 organs can be specified using `all`. |
| input_vcf          | The path to the input variants file |

### Config file

For details on setting up the config file for your machine, please read `config/README.md` and follow the instructions before continuing.

For your own runs, please make a copy of `config/config.yml` to your project directory and ensure that `runs_table` and `base_dir` are properly updated before running. 

### Running the example workflow

To run the workflow from the command line, set the working directory to:

```bash
cd path/to/TLand-predict
```

Set up the machine-specific configurations for `example/config.yml`.

To run the workflow with example files using **conda**:

```bash
snakemake -F --cores 1 --resources mem_mb=5000 gpu=1 queries=1 --configfile example/config.yml --use-conda
```

`--gpu` specifies the number of available GPUs to use to generate Sei variant effect prediction features. Can omit to use CPU only. 

`--queries` specifies the number of RegulomeDB queries that can be run in parallel. Values greater than 1 only work if multiple runs are submitted via the runs table. A query is performed for each TLand-predict run and can be quite memory-intensive if the variant lists are large. If scoring millions of variants, we recommend splitting the input variants into runs of 1 million variants each and set `--queries=1`. For multi-run workflows with smaller variant lists, feel free to bump up the number of concurrent queries to speed up the workflow.

For your own runs, the same command can be used. Just change `--configfile example/config.yml` to `--configfile /path/to/your/config.yml`. `-F` forces Snakemake to run from the beginning, so remove this option if you want to continue a run mid-way.

## References

> Zhao, N., Dong, S. & Boyle, A. P. Organ-specific prioritization and annotation of non-coding regulatory variants in the human genome. 2023.09.07.556700 Preprint at https://doi.org/10.1101/2023.09.07.556700 (2023).