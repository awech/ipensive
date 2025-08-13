"""
Created on 30-Mar-2022
Modified on 3-Jul-2025
@author: awech
"""

import logging
import argparse
import pandas as pd
from obspy import UTCDateTime as utc
from ipensive.array_processing import process_array
from ipensive import ipensive_utils as utils


def parse_args():
    """
    Parse command-line arguments for the script.
    
    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        epilog="e.g.: python back_populate.py -c <filename.yml> -s 202201010000 -e 202201020000"
    )

    parser.add_argument(
        "-c",
        "--config",
        type=str,
        help="Path to the config file (yml). If not specified, the default config will be used.",
        default=utils.get_config_file(),
    )
    parser.add_argument(
        "-s",
        "--starttime",
        type=str,
        help="Start time in UTC: YYYYMMDDHHMM",
        required=True,
    )
    parser.add_argument(
        "-e",
        "--endtime",
        type=str,
        help="End time in UTC: YYYYMMDDHHMM",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files. Default is False",
    )
    parser.add_argument(
        "-np",
        "--no-plot",
        action="store_true",
        help="Don't plot results. (e.g., if you just want to write out the CSV files)",
    )
    parser.add_argument(
        "-a",
        "--arrays",
        type=str,
        help="Comma-separated list of arrays to process (use _ for arrays with spaces).",
    )
    parser.add_argument(
        "-l", "--log", help="Specify log output file"
    )

    return parser.parse_args()


def run_backpopulate(config, T1, T2, OVERWRITE, ARRAYS):
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

    args = parse_args()  # Parse command-line arguments
    config_file = args.config
    config = utils.load_config(config_file)
    config["plot"] = not args.no_plot
    print(args)
    if args.arrays is not None:
        ARRAYS = list(args.arrays.replace("_"," ").split(","))
    else:
        ARRAYS = []

    utils.setup_logging(utc.utcnow(), config, arg_opt=args.log)
    my_log = logging.getLogger(__name__)

    run_backpopulate(config, args.starttime, args.endtime, args.overwrite, ARRAYS)
