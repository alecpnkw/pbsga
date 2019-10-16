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
        #iterate
        index_values = zip(df["N7_Index"],df["S5_Index"])
        index_dicts = [{"Index_1": x, "Index_2": y} for (x,y) in index_values]
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
