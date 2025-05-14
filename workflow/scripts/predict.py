from pathlib import Path
import pickle
import time
import argparse
from collections import defaultdict

import pandas as pd

def predict(organ, generic_features, organsp_features, gt_100_tf_chip, model_ls):
    
    pred_dict = defaultdict(defaultdict(list).copy)
    df_all = pd.read_parquet(generic_features)

    organsp_df = pd.read_parquet(organsp_features)

    df_all = pd.merge(df_all, 
                organsp_df.drop_duplicates(subset=['chrom', 'end']), 
                on=['chrom', 'end'], how='left')

    # make prediction
    if organ in gt_100_tf_chip:
        pred_dict[organ][model_ls[0][0]].append(model_ls[0][1].predict_proba(df_all)[:, 1]) # TLand
    else:
        pred_dict[organ][model_ls[1][0]].append(model_ls[1][1].predict_proba(df_all)[:, 1]) # TLand lightest

    return df_all, pred_dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--generic_features", help="Path to generic features file", type=str) 
    parser.add_argument("--organsp_features", help="Path to organ-specific features file", type=str)
    parser.add_argument("--organ", help="Organ name", type=str)
    parser.add_argument("--organ_list", help="List of organs TLand can make predictions for", type=str)
    parser.add_argument("--models_path", help="Path to directory containing model files", type=str)
    parser.add_argument("--out", help="Output file path", type=str)
    args = parser.parse_args()

    generic_features = Path(args.generic_features)
    organsp_features = args.organsp_features
    organ = args.organ
    organ_list = args.organ_list
    models_path = Path(args.models_path)
    outfile = args.out

    organ_ls = []
    with Path(organ_list).open() as file:
        for line in file:
            line = line.strip()
            organ_ls.append(line)
    allowed_organ_args = organ_ls.copy()
    allowed_organ_args.append("all")
    
    if args.organ not in allowed_organ_args:
        print(f"Error: Argument {args.organ} is not valid. Must be one of {allowed_organ_args}.")
        exit(1)

    tland_path = models_path / 'TLand_organSp.pickle'
    tland_lightest_path = models_path / 'TLand_organSp_lightest.pickle'
    model_ls = [['TLand', pickle.load(tland_path.open(mode='rb'))], 
            ['TLand_lightest', pickle.load(tland_lightest_path.open(mode='rb'))]]
    

    # Organs with greater than 100 TF-ChIP experiments (Use TLand for these, TLand-lightest for everything else)
    gt_100_tf_chip = ["epithelium",
                        "blood",
                        "bodily fluid",
                        "exocrine gland",
                        "endocrine gland",
                        "liver",
                        "lung",
                        "kidney",
                        "mammary gland",
                        "brain",
                        "connective tissue",
                        "skin of body",
                        "uterus"]

    start = time.time()
    print(f'Predicting TLand scores for {organ}...')
    
    df_all, pred_dict = predict(organ, generic_features, organsp_features, gt_100_tf_chip, model_ls)
    output_df = df_all[['chrom', 'end', 'ref', 'alt']].copy()
    if organ in gt_100_tf_chip:
        output_df[organ+"_"+model_ls[0][0]] = pred_dict[organ][model_ls[0][0]][0]
    else:
        output_df[organ+"_"+model_ls[1][0]] = pred_dict[organ][model_ls[1][0]][0]

    output_df.rename(columns={'end': 'pos'}, inplace=True)
    output_df.to_csv(outfile, sep='\t', index=None, compression='gzip')

    elapsed = time.time() - start
    print(f'Total time: {elapsed:.2f} seconds')
