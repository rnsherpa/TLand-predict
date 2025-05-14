from pathlib import Path
import json
import argparse

import pandas as pd
import numpy as np
import pyBigWig

def get_generic_features(input_json):
    '''Function to output RegDB generic features from a json object
    Args:
        input_json(json): a json object contains RegDB features
    Returns:
        out(list of str): RegDB features ('chrom','end','CHIP','Chromatin_accessibility','PWM','FOOTPRINT','QTL','PWM_matched','FOOTPRINT_matched')
    '''
    out = [input_json['chrom'],input_json['end']]
    feature_names = ['ChIP','Chromatin_accessibility','PWM','Footprint','QTL','PWM_matched','Footprint_matched',  'IC_matched_max', 'IC_max']
    out_features = [int(input_json['features'][feature]) for feature in feature_names] #convert to binary values
    out.extend(out_features)
    return out

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract generic features for variants.")
    parser.add_argument('--input_vcf', type=str, required=True, help='Path to input VCF file')
    parser.add_argument('--input_jsonl', type=str, required=True, help='Path to RegDB features JSONL file')
    parser.add_argument('--input_sei', type=str, required=True, help='Path to SEI features file')
    parser.add_argument('--dnase_sig_path', type=str, required=True, help='Path to quantile-normalized DNase signal bigWig files')
    parser.add_argument('--chip_sig_path', type=str, required=True, help='Path to quantile-normalized ChIP signal bigWig files')
    parser.add_argument('--out', type=str, required=True, help='Output file path')
    args = parser.parse_args()

    input_vcf = args.input_vcf
    input_jsonl = args.input_jsonl
    input_sei = args.input_sei
    dnase_sig_path = Path(args.dnase_sig_path)
    chip_sig_path = Path(args.chip_sig_path)
    outfile = args.out

    ### READ IN VARIANTS
    df_all = pd.read_csv(input_vcf, sep='\t', usecols=[0,1,2,3,4], names=['chrom','end','id','ref','alt'])
        
    ### GET GENERIC REGDB QUERY FEATURES
    RegDB_generic_features = []
    with open(input_jsonl) as f:
        for line in f:
            if line != '\n': # json.loads isn't happy if it encounters an empty line
                var_json = json.loads(line.strip())
                RegDB_generic_features.append(get_generic_features(var_json))
    # Chromatin accessibility renamed to DNASE bc model only recognizes DNASE
    df_all = pd.merge(df_all, pd.DataFrame(RegDB_generic_features,columns = ['chrom','end','CHIP','DNASE','PWM','FOOTPRINT','EQTL_2','PWM_matched','FOOTPRINT_matched','IC_matched_max', 'IC_max']).drop_duplicates(), 
                    left_on=['chrom', 'end'], right_on=['chrom', 'end'], how='left')
    
    ### GENERIC DNASE SIGNALS
    for DNase_sig_feature in ['DNase_var','DNase_quantile95','DNase_quantile1','DNase_quantile2','DNase_quantile3']:
        bw = pyBigWig.open(str(dnase_sig_path / f'{DNase_sig_feature}.bw'))
        df_all[DNase_sig_feature] = [np.round(np.array(bw.values(chrom,end-1,end))[0],4) for chrom,end in zip(df_all['chrom'],df_all['end'])] 
        bw.close()

    ### GENERIC CHIP SIGNALS
    for ChIP_sig_feature in ['ChIP_var','ChIP_quantile95','ChIP_quantile1','ChIP_quantile2','ChIP_quantile3']:
        bw = pyBigWig.open(str(chip_sig_path / f'{ChIP_sig_feature}.bw'))
        df_all[ChIP_sig_feature] = [np.round(np.array(bw.values(chrom,end-1,end))[0],4) for chrom,end in zip(df_all['chrom'],df_all['end'])] 
        bw.close()
    
    ### SEI SEQUENCE CLASSES
    sei_features = pd.read_csv(input_sei, sep='\t')
    sei_features.rename(columns={'pos':'end'}, inplace=True)
    sei_features.drop(['seqclass_max_absdiff', 'strand', 'id'], axis=1, inplace=True)

    df_all = df_all.merge(sei_features, how='left', left_on=['chrom', 'end', 'ref', 'alt'], 
                    right_on=['chrom', 'end', 'ref', 'alt'])
    df_all.iloc[:, -40:] = df_all.iloc[:, -40:].fillna(df_all.iloc[:, -40:].mean())
    df_all['max_abs_diff'] = df_all.iloc[:, -40:].abs().max(axis=1)

    # Save as parquet
    df_all.to_parquet(outfile, index=False)