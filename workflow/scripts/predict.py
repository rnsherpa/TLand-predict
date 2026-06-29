from pathlib import Path
import json
import pickle
import time
import argparse
from collections import defaultdict

import pandas as pd

def predict(organ, generic_features, organsp_features, selected_model_ls):
    
    pred_dict = defaultdict(defaultdict(list).copy)
    df_all = pd.read_parquet(generic_features)

    organsp_df = pd.read_parquet(organsp_features)

    df_all = pd.merge(df_all, 
                organsp_df.drop_duplicates(subset=['chrom', 'end']), 
                on=['chrom', 'end'], how='left')

    # make prediction
    for model_name, model in selected_model_ls:
        pred_dict[organ][model_name].append(model.predict_proba(df_all)[:, 1])

    return df_all, pred_dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--generic_features", help="Path to generic features file", type=str) 
    parser.add_argument("--organsp_features", help="Path to organ-specific features file", type=str)
    parser.add_argument("--organ", help="Organ name", type=str)
    parser.add_argument('--organ_mapping_json', type=str, required=True, help='Path to file mapping organ names with underscores to organ names with spaces')
    parser.add_argument("--organ_list", help="List of organs TLand can make predictions for", type=str)
    parser.add_argument("--models_path", help="Path to directory containing model files", type=str)
    parser.add_argument("--all_TLand_models", help="Use all TLand models for prediction", action="store_true")
    parser.add_argument("--out", help="Output file path", type=str)
    args = parser.parse_args()

    generic_features = Path(args.generic_features)
    organsp_features = args.organsp_features
    safe_organ_name = args.organ
    organ_list = args.organ_list
    models_path = Path(args.models_path)
    all_TLand_models = args.all_TLand_models
    outfile = args.out

    # Get organ name with spaces since that is how they are labeled in RegDB JSON and bigwig filenames
    with open(args.organ_mapping_json, "r") as f:
        organ_mapping = json.load(f)
    organ = organ_mapping[safe_organ_name]

    organ_ls = []
    with Path(organ_list).open() as file:
        for line in file:
            line = line.strip()
            organ_ls.append(line)
    allowed_organ_args = organ_ls.copy()
    allowed_organ_args.append("all")
    
    if organ not in allowed_organ_args:
        print(f"Error: Argument {organ} is not valid. Must be one of {allowed_organ_args}.")
        exit(1)

    tland_path = models_path / 'TLand_organSp.pickle'
    tland_light_path = models_path / 'TLand_organSp_light.pickle'
    tland_lightest_path = models_path / 'TLand_organSp_lightest.pickle'
    model_ls = [['TLand', pickle.load(tland_path.open(mode='rb'))], 
            ['TLand_light', pickle.load(tland_light_path.open(mode='rb'))],
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
    
    if all_TLand_models:
        selected_model_ls = model_ls
    elif organ in gt_100_tf_chip:
        selected_model_ls = [model_ls[0]]
    else:
        selected_model_ls = [model_ls[1]]

    df_all, pred_dict = predict(organ, generic_features, organsp_features, selected_model_ls)
    output_df = df_all[['chrom', 'end', 'id', 'ref', 'alt']].copy()
    for model_name, _ in selected_model_ls:
        output_df[safe_organ_name+"."+model_name] = pred_dict[organ][model_name][0]

    output_df.rename(columns={'end': 'pos'}, inplace=True)
    output_df.to_csv(outfile, sep='\t', index=None, compression='gzip')

    elapsed = time.time() - start
    print(f'Total time: {elapsed:.2f} seconds')
