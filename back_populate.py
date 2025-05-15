#!/home/rtem/miniconda3/envs/ipensive/bin/python
# -*- coding: utf-8 -*-
"""
Created on 30-Mar-2022
@author: awech
"""

from pathlib import Path
import pandas as pd
from obspy import UTCDateTime
from array_processing import process_array, parse_args
import ipensive_utils as utils


T1 = "2025-05-15 15:40"
T2 = "2025-05-15 18:40"
OVERWRITE = True
PLOT = True

# specify ARRAYS if you want to just process specific arrays
ARRAYS = []
# ARRAYS = ["Wake Island North", "Wake Island South"]

args = parse_args()
config = utils.load_config(args.config)
config["plot"] = PLOT

def get_file_path(t, array_name, config):

    out_web_dir = Path(config["OUT_WEB_DIR"])
    array_dict = config.get(array_name, {})
    network_dir = out_web_dir / array_dict["NETWORK_NAME"]
    array_dir = network_dir / array_dict["ARRAY_NAME"]
    year_dir = array_dir / str(t.year)
    julian_day_dir = year_dir / str(t.day_of_year)

    time = t.strftime("%Y%m%d-%H%M")
    file = julian_day_dir / f"{array_name}_{time}.png"

    return file


def run_backpopulate():
    t1 = UTCDateTime(T1) + config["PARAMS"]["DURATION"]
    for t in pd.date_range(T2, t1.strftime("%Y%m%d%H%M"), freq="-10min"):
        print(t)

        if len(ARRAYS) > 0:
            array_list = ARRAYS
        else:
            array_list = config["array_list"]
        for array_name in array_list:
            # check if you should process this time window
            file = get_file_path(t, array_name, config)
            if not file.exists() or OVERWRITE:
                process_array(config, array_name, UTCDateTime(t))
                print("\n")
            else:
                print("File exists. No overwrite. Skip " + array_name)
    return


if __name__ == '__main__':
    run_backpopulate()
