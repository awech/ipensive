"""
Created on 30-Mar-2022
Modified on 10-Jun-2025
@author: awech
"""

import logging
import pandas as pd
from obspy import UTCDateTime as utc
from ipensive.array_processing import process_array, parse_args
from ipensive import ipensive_utils as utils


T1 = "2025-06-10 10:00"
T2 = "2025-06-10 11:00"
OVERWRITE = False
PLOT = True

# specify ARRAYS if you want to just process specific arrays
ARRAYS = []
# ARRAYS = ["Wake Island North", "Wake Island South"]

args = parse_args()
config = utils.load_config(args.config)
config["plot"] = PLOT

utils.setup_logging(utc.utcnow(), config, arg_opt=args.log)
my_log = logging.getLogger(__name__)

def run_backpopulate():
    t1 = utc(T1) + config["PARAMS"]["DURATION"]
    for t in pd.date_range(T2, t1.strftime("%Y%m%d%H%M"), freq="-10min"):
        my_log.info(t)

        if len(ARRAYS) > 0:
            array_list = ARRAYS
        else:
            array_list = config["array_list"]
        for array_name in array_list:
            # check if you should process this time window
            file = utils.get_file_path(t, array_name, config)
            if not file.exists() or OVERWRITE:
                process_array(config, array_name, utc(t))
                my_log.info("\n")
            else:
                my_log.info("File exists. No overwrite. Skip " + array_name)


if __name__ == '__main__':
    run_backpopulate()
