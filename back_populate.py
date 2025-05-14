#!/home/rtem/miniconda3/envs/ipensive/bin/python
# -*- coding: utf-8 -*-
"""
Created on 30-Mar-2022
@author: awech
"""

import os
import pandas as pd
from obspy import UTCDateTime
from array_processing import process_array, parse_args
import ipensive_utils as utils

T1 = "2025-05-12 00:00"
T2 = "2025-05-14 21:50"
OVERWRITE = True

# specify ARRAYS if you want to just process specific arrays
ARRAYS = []
# ARRAYS = ["Wake Island North", "Wake Island South"]

args = parse_args()
config = utils.load_config(args.config)

def make_file_path(t, array_name, config):
    year = "{:.0f}".format(t.year)
    day = "{:03.0f}".format(t.dayofyear)
    time = t.strftime("%Y%m%d-%H%M")
    array_dict = config[array_name]
    file = "{}/{}/{}/{}/{}/{}_{}.png".format(
        config["OUT_WEB_DIR"],
        array_dict["NETWORK_NAME"],
        array_dict["ARRAY_NAME"],
        year,
        day,
        array_dict["ARRAY_NAME"],
        time,
    )
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
            file = make_file_path(t, array_name, config)
            if not os.path.exists(file) or OVERWRITE:
                process_array(config, array_name, UTCDateTime(t))
                print("\n")
            else:
                print("File exists. No overwrite. Skip " + array_name)
    return


if __name__ == '__main__':
    run_backpopulate()
