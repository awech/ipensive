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
        required=False,
    )
    parser.add_argument(
        "-e",
        "--endtime",
        type=str,
        help="End time in UTC: YYYYMMDDHHMM",
        required=False,
    )
    parser.add_argument(
        "-dt",
        "--duration",
        type=str,
        help="Duration in hours (\"h\") days (\"d\") before present to process (e.g. -dt 3h)",
        required=False,
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
    """
    Run the back population process for the specified time window and arrays.

    Args:
        config (dict): Configuration dictionary.
        T1 (str): Start time in UTC.
        T2 (str): End time in UTC.
        OVERWRITE (bool): Flag to overwrite existing files.
        ARRAYS (list): List of arrays to process.
    """

    t1 = utc(T1) + config["PARAMS"]["DURATION"]
    for t in pd.date_range(T2, t1.strftime("%Y%m%d%H%M"), freq="-10min"):
        my_log.info(t)

        if len(ARRAYS) > 0:
            array_list = ARRAYS
        else:
            array_list = config["array_list"]
        for array_name in array_list:
            # check if you should process this time window
            file = utils.get_pngfile_path(t, array_name, config)
            if not file.exists() or OVERWRITE:
                try:
                    process_array(config, array_name, utc(t))
                except Exception as ex:
                    my_log.error(f"Error processing {array_name}: {ex}")
                my_log.info("\n")
            else:
                my_log.info("File exists. No overwrite. Skip " + array_name)
    my_log.info("Back population complete.")
    my_log.info("Writing .html file")
    utils.write_html(config)


if __name__ == '__main__':
    """
    Main entry point for the back population script.
    """

    args = parse_args()  # Parse command-line arguments
    ipensive_config_file = args.config
    config = utils.load_ipensive_config(ipensive_config_file)
    config["plot"] = not args.no_plot

    if args.arrays is not None:
        ARRAYS = list(args.arrays.replace("_"," ").split(","))
    else:
        ARRAYS = []

    utils.setup_logging(utc.utcnow(), config, arg_opt=args.log)
    my_log = logging.getLogger(__name__)
    for f in config["array_config_files"]:
        my_log.info(f"Array config: {f}")
    for key, value in args.__dict__.items():
        if key in ["starttime", "endtime"]:
            if value is not None:
                value = utc(value).strftime("%Y-%m-%d %H:%M:%S")
        my_log.info(f"{key}: {value}")
    my_log.info("\n")

    if args.starttime is not None and args.endtime is not None:
        T1 = args.starttime
        T2 = args.endtime
        if utc(T1) > utc(T2):
            raise ValueError("Start time must be before end time.")
    elif args.duration is not None:
        duration_value = args.duration[:-1]
        if args.duration.endswith("h"):
            dt = pd.Timedelta(hours=int(duration_value))
        elif args.duration.endswith("d"):
            dt = pd.Timedelta(days=int(duration_value))
        else:
            raise ValueError("Invalid duration format. Use 'h' for hours or 'd' for days.")

        if args.starttime is not None:
            T1 = utc(args.starttime)
            T2 = T1 + dt
        elif args.endtime is not None:
            T2 = utc(args.endtime)
            T1 = T2 - dt
        else:
            T2 = utc.utcnow()
            T1 = T2 - dt
            my_log.info(f"No start or end time specified. Using current time {T2.strftime('%Y-%m-%d %H:%M:%S')} as end time")

    if "T1" not in locals() and "T2" not in locals():
        raise ValueError("Must define start (-s) and end (-e) times, or one start/end with a duration (-dt)")

    date_fmt = "%Y%m%d%H%M"
    T1 = utc(T1).strftime(date_fmt)[:-1] + "0"
    T2 = utc(T2).strftime(date_fmt)[:-1] + "0"

    run_backpopulate(config, T1, T2, args.overwrite, ARRAYS)
