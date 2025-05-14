from pathlib import Path
import argparse
import json
from collections import defaultdict

import numpy as np
import pandas as pd
import pyBigWig

def get_organ_sp_RegDB_features(var_json, organ):
    '''Function to output RegDB organ specific features from a json object of a variant
    Args:
        var_json(json): a json object contains RegDB features & peaks information
        organ (str): biosample organ name (e.g. ['bodily fluid','blood','brain'])
    Returns:
        feature_values (list): values for 'chrom','end','DNASE','FOOTPRINT', 'ChIP-seq', 'unique targets', 'unique biosamples ChIP', 'unique biosamples DNase'
    '''
    method_names = ['DNase-seq', 'footprints', 'ChIP-seq']
    results = defaultdict(int) #keys: method_names; default value is 0
    unique_targets = []
    unique_biosamples = []
    unique_biosamples_dnase = []
    for exp in var_json['peaks']:
        try:    
            if exp['method'] in method_names:
                try:
                    if exp['method'] == 'ChIP-seq' and exp['targets'][0].startswith('POLR'): #not counting POLR ChIP-seq
                        continue
                except IndexError: #some json does not have an empty 'targets', and we still count it
                    exp['targets'] = ['NA']
                if organ in exp['organ_slims']:
                    results[exp['method']] += 1
                    if exp['method'] == 'ChIP-seq':
                        unique_targets.append(exp['targets'][0])
                        unique_biosamples.append(exp['biosample_term_name'])
                    if exp['method'] == 'DNase-seq':
                        unique_biosamples_dnase.append(exp['biosample_term_name'])
        except KeyError:
            continue
    feature_values = [var_json['chrom'],var_json['end']]
    feature_values.append(results['DNase-seq'])
    feature_values.append(results['footprints'])
    feature_values.append(results['ChIP-seq'])
    feature_values.append(len(set(unique_targets))) #unique count TFs in ChIP
    feature_values.append(len(set(unique_biosamples))) #unique count biosamples in ChIP
    feature_values.append(len(set(unique_biosamples_dnase))) #unique count biosamples in DNase
    return feature_values

def get_organ_sp_RegDB_features_ctcf(var_json, organ):
    '''Function to output RegDB organ specific ctcf features from a json object of a variant
    Args:
        var_json(json): a json object contains RegDB features & peaks information
        organ (str): biosample organ name (e.g. ['bodily fluid','blood','brain'])
    Returns:
        feature_values (list): values for 'chrom','end','CTCF','unique biosamples ChIP'
    '''
    method_names = ['ChIP-seq']
    results = defaultdict(int) #keys: method_names; default value is 0
    unique_biosamples = []
    for exp in var_json['peaks']:
        try:
            if exp['method'] in method_names:
                try:
                    if exp['method'] == 'ChIP-seq' and exp['targets'][0].startswith('CTCF'):
                        if organ in exp['organ_slims']:
                            results['CTCF'] += 1
                            unique_biosamples.append(exp['biosample_term_name'])
                except IndexError:
                    continue
        except KeyError:
            continue
    feature_values = [var_json['chrom'],var_json['end']]
    feature_values.append(results['CTCF'])
    feature_values.append(len(set(unique_biosamples))) #unique count biosamples in CTCF
    return feature_values

def get_organ_sp_features_histone(var_json, organ):
    histone_list = ['H3K27ac','H3K36me3','H3K4me1','H3K4me3','H3K27me3']
    results = defaultdict(int) #keys: histone mark; default value is 0
    for exp in var_json['peaks']:
        try:    
            if exp['method'] == 'Histone ChIP-seq':
                if organ in exp['organ_slims']:
                    results[exp['target_label']] += 1
        except KeyError:
            continue
    feature_values = [var_json['chrom'],var_json['end']]
    for histone in histone_list:
        feature_values.append(results[histone])
    return feature_values

if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Extract organ-specific features for variants.")
    parser.add_argument('--input_jsonl', type=str, required=True, help='Path to RegDB features JSONL file')
    parser.add_argument("--organ", help="Organ name", type=str) 
    parser.add_argument('--organsp_dnase_sig_path', type=str, required=True, help='Path to quantile-normalized organ-specific DNase signal bigWig files')
    parser.add_argument('--total_num_path', type=str, required=True, help='Path to files containing total number of organ-specific annotations for a given feature')
    parser.add_argument('--out', type=str, required=True, help='Output file path')
    args = parser.parse_args()

    input_jsonl = args.input_jsonl
    organ = args.organ
    organsp_dnase_sig_path = Path(args.organsp_dnase_sig_path)
    total_num_path = Path(args.total_num_path)
    outfile = args.out

    total_num_dict={}
    for feature in ['DNASE', 'TF', 'CTCF', 'H3K27ac', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K27me3']:
        total_num_dict[feature] = dict(pd.read_csv(total_num_path / f'{feature}_totalNum_organ_hg38.txt', sep='\t', header=None).values)

    pseudo_count = 2

    ### DNase, footprint, and ChIP features
    regDB_organSp_features = []
    with open(input_jsonl) as f:
        print(f'Getting DNase, footprint, and ChIP feature counts for {organ} from RegDB query...')
        for line in f:
            if line != '\n':
                var_json = json.loads(line.strip())
                var_features = get_organ_sp_RegDB_features(var_json, organ)
                regDB_organSp_features.append(var_features)

    organSp_df = pd.DataFrame(regDB_organSp_features, columns = [
        'chrom','end', 'DNASE_organSp', 'FOOTPRINT_organSp', 'CHIP_organSp', 
        'CHIP_organSp_uniq', 'CHIP_organSp_biosample_uniq', 'DNASE_organSp_biosample_uniq'])
    try:
        organSp_df['CHIP_organSp_perc'] = organSp_df['CHIP_organSp'] / total_num_dict['TF'][organ]
    except KeyError:
        organSp_df['CHIP_organSp_perc'] = organSp_df['CHIP_organSp'] / pseudo_count
    try:
        organSp_df['DNASE_organSp_perc'] = organSp_df['DNASE_organSp'] / total_num_dict['DNASE'][organ]
    except KeyError:
        organSp_df['DNASE_organSp_perc'] = organSp_df['DNASE_organSp'] / pseudo_count

    ### CTCF ChIP features
    regDB_organSp_features_ctcf = []
    with open(input_jsonl) as f:
        print(f'Getting CTCF feature counts for {organ} from RegDB query...')
        for line in f:
            if line != '\n':
                var_json = json.loads(line.strip())
                var_features = get_organ_sp_RegDB_features_ctcf(var_json, organ)
                regDB_organSp_features_ctcf.append(var_features)

    organSp_df = organSp_df.merge(pd.DataFrame(regDB_organSp_features_ctcf, columns = ['chrom','end','CTCF_organSp','CTCF_organSp_biosample_uniq']), on=['chrom', 'end'], how='inner')
    
    try:
        organSp_df['CTCF_organSp_perc'] = organSp_df['CTCF_organSp'] / total_num_dict['CTCF'][organ]
    except:
        organSp_df['CTCF_organSp_perc'] = organSp_df['CTCF_organSp'] / pseudo_count

    ### Histone ChIP features
    regDB_organSp_features_histone = []
    with open(input_jsonl) as f:
        print(f'Getting histone feature counts for {organ} from RegDB query...')
        for line in f:
            if line != '\n':
                var_json = json.loads(line.strip())
                var_features = get_organ_sp_features_histone(var_json, organ)
                regDB_organSp_features_histone.append(var_features)
        
    organSp_df = organSp_df.merge(pd.DataFrame(regDB_organSp_features_histone, columns = [
        'chrom','end','H3K27ac_organSp','H3K36me3_organSp','H3K4me1_organSp','H3K4me3_organSp','H3K27me3_organSp']))
    
    try:
        organSp_df['H3K27ac_organSp_perc'] = organSp_df['H3K27ac_organSp'] / total_num_dict['H3K27ac'][organ]
    except KeyError:
        organSp_df['H3K27ac_organSp_perc'] = organSp_df['H3K27ac_organSp'] / pseudo_count
    
    try:
        organSp_df['H3K36me3_organSp_perc'] = organSp_df['H3K36me3_organSp'] / total_num_dict['H3K36me3'][organ]
    except KeyError:
        organSp_df['H3K36me3_organSp_perc'] = organSp_df['H3K36me3_organSp'] / pseudo_count
    
    try:
        organSp_df['H3K4me1_organSp_perc'] = organSp_df['H3K4me1_organSp'] / total_num_dict['H3K4me1'][organ]
    except KeyError:
        organSp_df['H3K4me1_organSp_perc'] = organSp_df['H3K4me1_organSp'] / pseudo_count    
    
    try:
        organSp_df['H3K4me3_organSp_perc'] = organSp_df['H3K4me3_organSp'] / total_num_dict['H3K4me3'][organ]
    except KeyError:
        organSp_df['H3K4me3_organSp_perc'] = organSp_df['H3K4me3_organSp'] / pseudo_count
    
    try:
        organSp_df['H3K27me3_organSp_perc'] = organSp_df['H3K27me3_organSp'] / total_num_dict['H3K27me3'][organ]
    except KeyError:
        organSp_df['H3K27me3_organSp_perc'] = organSp_df['H3K27me3_organSp'] / pseudo_count
    
    ### Signal features
    for DNase_sig_feature in ['DNase_var','DNase_quantile95','DNase_quantile1','DNase_quantile2','DNase_quantile3']:
        bw = pyBigWig.open(str(organsp_dnase_sig_path / f'{DNase_sig_feature}_{organ}.bw'))
        organSp_df[f'{DNase_sig_feature}_organSp'] = [np.round(np.array(bw.values(chrom,end-1,end))[0],4) for chrom,end in zip(organSp_df['chrom'],organSp_df['end'])] 
        bw.close()

    organSp_df.drop_duplicates(subset=['chrom', 'end']).to_parquet(outfile, index=None)