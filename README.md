# TLand-predict

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)

A Snakemake pipeline to create input feature tables and run predictions with TLand

## Usage

Detailed information about input data and workflow configuration can also be found in the [`config/README.md`](config/README.md).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository or its DOI.

## Deployment options

To run the workflow from command line, change the working directory.

```bash
cd path/to/TLand-predict
```

Adjust options in the default config file `config/config.yml`.
Before running the complete workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

To run the workflow with test files using **conda**:

```bash
snakemake -F --cores 1 --resources mem_mb=1000 --configfile tests/config.yml --use-conda
```

## References

> Zhao, N., Dong, S. & Boyle, A. P. Organ-specific prioritization and annotation of non-coding regulatory variants in the human genome. 2023.09.07.556700 Preprint at https://doi.org/10.1101/2023.09.07.556700 (2023).