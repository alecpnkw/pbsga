import pandas as pd
import yaml, sys, os

def create_config(input_csv_fhs, target_size = 2100, index_type = "nextera"):
    datasets = []
    configs = []
    for input_csv in input_csv_fhs:
        #read in primer sheet
        df = pd.read_csv(input_csv)
        #filter to defined N7/S5
        df = df[(pd.notna(df["N7_Index"])) & (pd.notna(df["S5_Index"]))]
        sub_df = df[["N7_Index","S5_Index","Forward_Primer_2ndRd_Sequence","Reverse_Primer_2ndRd_Sequence"]]
        sub_df = sub_df.assign(fwd_index = sub_df["N7_Index"],
                     rev_index = sub_df["S5_Index"],
                     fwd_primer = sub_df["Forward_Primer_2ndRd_Sequence"],
                     rev_primer = sub_df["Reverse_Primer_2ndRd_Sequence"]
                     )
        sub_df = sub_df[["fwd_index","rev_index","fwd_primer","rev_primer"]]
        #iterate
        index_dicts = [sub_df.iloc[i].to_dict() for i in range(sub_df.shape[0])]
        template_dict = dict(zip(df["Sample"], index_dicts))
        dataset_dict = {
            "target_size": target_size,
            "index_type": index_type,
            "templates": template_dict
            }
        configs.append(dataset_dict)
        datasets.append(os.path.basename(input_csv).split(".")[0])

    sys.stdout.write(yaml.dump({"datasets": dict(zip(datasets, configs))}))

create_config(sys.argv[1:])
