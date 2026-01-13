# Setting up the config

Snakemake uses config files to specify things like input files and workflow parameters. This document goes over the necessary steps to set up the config file for `TLand-predict`.

### Machine-specific configurations

Running `TLand-predict` requires prior setup of specific software and data to your machine. This will result in paths unique to your setup being added to the `config.yml` file. These are a ONE TIME edit to `config.yml` and should remain unchanged for any runs performed on your machine.

1. Set up a local RegulomeDB server following the instructions from [this repository](https://github.com/ENCODE-DCC/genomic-data-service)
- Set the `gds_dir` key in the config to the path to `genomic-data-service`

2. Download the necessary Sei files from:
- Set the `sei_dir` key in the config to the path to `Sei`

3. Download the remaining required files from [] and set:
- `dnase_sig_path` to 
- `chip_sig_path` to
- `organsp_dnase_sig_path` to
- `models_path` to

Configs under `# Resources` do not require changing as the files are in the repository.

### Run-specific configurations

`runs_table`: Path to runs table

`base_dir`: Path to directory where runs should be output