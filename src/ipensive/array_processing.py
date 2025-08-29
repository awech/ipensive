"""
Created on 24-Apr-2018
Modified on 15-May-2025
@author: awech
"""

import os
import argparse
import logging
import time
import pandas as pd
from obspy import UTCDateTime as utc
from . import ipensive_utils as utils
from .plotting_utils import plot_results, save_figure
from . import data_utils
from .metadata_utils import add_metadata

from lts_array import ltsva
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, append=True)

my_log = logging.getLogger(__name__)

def parse_args():
    """
    Parse command-line arguments for the script.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        epilog="e.g.: python array_processing.py -c <filename.yml>"
    )

    parser.add_argument(
        "-c",
        "--config",
        type=str,
        help="Name of the config file (yml)",
        default=utils.get_config_file(),
    )
    parser.add_argument(
        "-t",
        "--time",
        type=str,
        help="UTC time stamp: YYYYMMDDHHMM (optional, otherwise grabs current UTC time)",
    )
    parser.add_argument(
        "-a",
        "--array",
        type=str,
        help="Array name if you want to process a single array. (Use _ instead of spaces, if necessary)",
    )
    parser.add_argument(
        "-np",
        "--no-plot",
        action="store_true",
        help="Don't plot results. (e.g., if you just want to write out the CSV files)",
    )
    parser.add_argument(
        "-l",
        "--log",
        type=str,
        help="Log file name (optional)",
    )

    return parser.parse_args()


def get_starttime(config, args):
    """
    Determine the start time for processing based on arguments or current time.

    Args:
        config (dict): Configuration dictionary.
        args (argparse.Namespace): Parsed command-line arguments.

    Returns:
        tuple: Start time (obspy.UTCDateTime) and delay (int).
    """
    date_fmt = "%Y-%m-%d %H:%M"

    if args.time:
        # Use the provided time argument
        T0 = utc(args.time)
        wait = 0
    else:
        # Use the current UTC time and add latency and window length
        T0 = utc.utcnow()
        wait = config["PARAMS"]["LATENCY"] + config["PARAMS"]["WINDOW_LENGTH"]

    # Round down to the nearest 10-minute interval
    T0 = utc(T0.strftime(date_fmt)[:-1] + "0")
    earliest_start = T0 + wait
    delay = max(0, earliest_start - utc.utcnow())

    return T0, delay


def do_LTS(st, array_params, lat_list, lon_list, skip_chans):

    my_log.info("Performing LTS analysis...")
    overlap_fraction = array_params["OVERLAP"] / array_params["WINDOW_LENGTH"]
    ALPHA = array_params["LTS_ALPHA"] if len(st) > 3 else 1.0
    skip_inds = [i for i, tr in enumerate(st) if tr.id in skip_chans]
    velocity, azimuth, t, mccm, lts_dict, sigma_tau, Vel_err, Baz_err = ltsva(
        st.copy(), lat_list, lon_list, array_params["WINDOW_LENGTH"], overlap_fraction, alpha=ALPHA, remove_elements=skip_inds
    )

    df = pd.DataFrame({
        "Time": t,
        "Array": array_params["ARRAY_NAME"],
        "Azimuth": azimuth,
        "Velocity": 1000 * velocity, # Convert velocity to m/s
        "MCCM": mccm,
        "Pressure": data_utils.get_pressures(st, t, array_params),
        "Sigma_tau": sigma_tau,
        "Vel_err": 1000 * Vel_err, # Convert to m/s
        "Baz_err": Baz_err
    })

    return df, lts_dict


def process_array(config, array_name, T0, return_figure=False):
    """
    Process a single array for the specified time window.

    Args:
        config (dict): Configuration dictionary.
        array_name (str): Name of the array to process.
        T0 (obspy.UTCDateTime): End time of the processing window.

    Returns:
        matplotlib figure or None
    """

    array_params = config[array_name]
    t1 = T0 - array_params["DURATION"]  # Start time of the processing window
    t2 = T0  # End time of the processing window

    T1 = t1 - array_params["TAPER"]  # Extended start time for tapering
    T2 = t2 + array_params["WINDOW_LENGTH"] + array_params["TAPER"]  # Extended end time

    my_log.info("--- " + array_params["ARRAY_NAME"] + " ---")
    my_log.info(f"{t1.strftime('%Y-%b-%d %H:%M')} - {t2.strftime('%H:%M')}")

    # Check for minimum channels
    if len(array_params["NSLC"]) < array_params["MIN_CHAN"]:
        my_log.warning("Not enough channels defined.")
        return None
    # Check for latency delay
    if os.getenv("FROMCRON") == "yep":
        if "EXTRA_PAUSE" in array_params:
            time.sleep(array_params["EXTRA_PAUSE"])  # Pause if running from a cron job

    # Download data
    st = data_utils.grab_data(array_params["CLIENT"], array_params["NSLC"], T1, T2)
    # st.write("test_raw.mseed")
    
    # Check data quality
    good_data, skip_chans = data_utils.QC_data(st, array_params)
    if not good_data:
        return None

    # Add metadata or coordinates
    st, lat_list, lon_list = add_metadata(st, config, array_name, skip_chans)

    # Add volcano backazimuths
    array_params = utils.get_target_backazimuth(st, config, array_params)

    # Preprocess data
    st = data_utils.preprocess_data(st, t1, t2, skip_chans, array_params)
    # st.write("test_preprocessed.mseed")

    # Perform array processing
    results_df, lts_dict = do_LTS(st, array_params, lat_list, lon_list, skip_chans)
    # results_df.to_csv("test_results.csv", index=False, float_format='%.7f')

    # Write output files
    utils.write_data_files(t2, st, results_df, config)

    # Generate plots
    if not config["plot"]:
        return None

    if not return_figure:
        utils.web_folders(t2, config, array_params)    
        for plotsize in ["big", "thumbnail"]:
            fig = plot_results(t2, st, results_df, lts_dict, skip_chans, array_params, plotsize)
            save_figure(fig, config, array_params, t2, plotsize)
    else:
        fig = plot_results(t2, st, results_df, lts_dict, skip_chans, array_params, "big")

    return fig