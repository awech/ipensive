#!/home/rtem/miniconda3/envs/ipensive/bin/python
# -*- coding: utf-8 -*-
"""
Created on 24-Apr-2018
Modified on 15-May-2025
@author: awech
"""

import os
import sys
sys.dont_write_bytecode = True  # Prevent Python from creating .pyc files (avoids clutter)
import numpy as np
import jinja2
from obspy.core import UTCDateTime as utc
from matplotlib import dates
import time
import ipensive_utils as utils
import argparse
from lts_array import ltsva
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, append=True)

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
        default="config.yml",
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
        "--no-plot",
        action="store_true",
        help="Disable plotting of results",
    )

    return parser.parse_args()


def get_starttime(config, args):
    """
    Determine the start time for processing based on arguments or current time.

    Args:
        config (dict): Configuration dictionary.
        args (argparse.Namespace): Parsed command-line arguments.

    Returns:
        tuple: Start time (UTCDateTime) and delay (int).
    """
    date_fmt = "%Y-%m-%d %H:%M"

    if args.time:
        # Use the provided time argument
        T0 = utc(args.time)
        delay = 0
    else:
        # Use the current UTC time and add latency and window length
        T0 = utc.utcnow()
        delay = config["PARAMS"]["LATENCY"] + config["PARAMS"]["WINDOW_LENGTH"]
        print(f"Waiting {delay:g} seconds")

    # Round down to the nearest 10-minute interval
    T0 = utc(T0.strftime(date_fmt)[:-1] + "0")

    return T0, delay


def write_html(config):
    """
    Generate an HTML file for web output.

    Args:
        config (dict): Configuration dictionary.

    Returns:
        None
    """
    if "EXTRA_LINKS" not in config.keys():
        config["EXTRA_LINKS"] = []
    script_path = os.path.dirname(__file__)
    with open(os.path.join(script_path, "index.template"), "r") as f:
        template = jinja2.Template(f.read())
    html = template.render(
        networks=config["network_list"],
        arrays=config["NETWORKS"],
        extra_links=config["EXTRA_LINKS"],
    )
    with open(os.path.join(config["OUT_WEB_DIR"], "index.html"), "w") as f:
        f.write(html)


def process_array(config, array_name, T0):
    """
    Process a single array for the specified time window.

    Args:
        config (dict): Configuration dictionary.
        array_name (str): Name of the array to process.
        T0 (UTCDateTime): End time of the processing window.

    Returns:
        None
    """
    array_params = config[array_name]
    t1 = T0 - array_params["DURATION"]  # Start time of the processing window
    t2 = T0  # End time of the processing window

    T1 = t1 - array_params["TAPER"]  # Extended start time for tapering
    T2 = t2 + array_params["WINDOW_LENGTH"] + array_params["TAPER"]  # Extended end time

    print("--- " + array_params["ARRAY_NAME"] + " ---")
    print(f"{t1.strftime('%Y-%b-%d %H:%M')} - {t2.strftime('%H:%M')}")
    if os.getenv("FROMCRON") == "yep":
        time.sleep(array_params["EXTRA_PAUSE"])  # Pause if running from a cron job

    #### Download data ####
    if len(array_params["NSLC"]) < array_params["MIN_CHAN"]:
        print("Not enough channels defined.")
        return

    st = utils.grab_data(
        array_params["CLIENT"],
        array_params["NSLC"],
        T1,
        T2,
    )

    #### Check for enough data ####
    check_st = st.copy()
    skip_chans = []
    for tr in check_st:
        if np.sum(np.abs(tr.data)) == 0:  # Check for blank traces
            skip_chans.append(tr.id)
            check_st.remove(tr)
    if len(check_st) < array_params["MIN_CHAN"]:
        print("Too many blank traces. Skipping.")
        return
    ########################

    #### Check for gappy data ####
    for tr in check_st:
        if np.any([np.any(tr.data == 0)]):  # Check for gaps in data
            check_st.remove(tr)
    if len(check_st) < array_params["MIN_CHAN"]:
        print("Too gappy. Skipping.")
        return
    ########################

    #### Add metadata or coordinates ####
    if isinstance(array_params["NSLC"], dict):
        st = utils.add_coordinate_info(st, config, array_name)
    else:
        st = utils.add_metadata(st, config, skip_chans)
    array_params = utils.get_target_backazimuth(st, config, array_params)
    ########################

    #### Preprocess data ####
    for tr in st:
        if tr.id in skip_chans:
            continue
        if isinstance(array_params["NSLC"], dict):
            tr.data = tr.data / array_params["NSLC"][tr.id]["gain"]
        else:
            tr.remove_sensitivity(tr.inventory)
    st.detrend("demean")
    for tr in st:
        if tr.stats["sampling_rate"] == 100.0:
            tr.decimate(2)
        if tr.stats["sampling_rate"] != 50.0:
            tr.resample(50.0)
        if tr.stats["sampling_rate"] == 50.0:
            if float(array_params["FREQMAX"]) < 50 / 4.0:
                tr.decimate(2)
    st.taper(max_percentage=None, max_length=array_params["TAPER"])
    st.filter(
        "bandpass",
        freqmin=array_params["FREQMIN"],
        freqmax=array_params["FREQMAX"],
        corners=2,
        zerophase=True,
    )
    st.trim(t1, t2 + array_params["WINDOW_LENGTH"])
    ########################

    #### Perform array processing ####
    lat_list = []
    lon_list = []
    for tr in st:
        lat_list.append(tr.stats.coordinates.latitude)
        lon_list.append(tr.stats.coordinates.longitude)
    overlap_fraction = array_params["OVERLAP"] / array_params["WINDOW_LENGTH"]
    ALPHA = array_params["LTS_ALPHA"] if len(st) > 3 else 1.0
    skip_inds = [i for i, tr in enumerate(st) if tr.id in skip_chans]
    velocity, azimuth, t, mccm, lts_dict, sigma_tau, *_ = ltsva(
        st.copy(), lat_list, lon_list, array_params["WINDOW_LENGTH"], overlap_fraction, alpha=ALPHA, remove_elements=skip_inds
    )
    pressure = []
    for tr_win in st[0].slide(
        window_length=array_params["WINDOW_LENGTH"],
        step=array_params["WINDOW_LENGTH"] - array_params["OVERLAP"],
    ):
        pressure.append(np.max(np.abs(tr_win.data)))
    pressure = np.array(pressure)

    #### Generate plots ####
    if config["plot"]:
        try:
            print("Setting up web output folders")
            utils.web_folders(t2, config, array_params)
            print("Making plot...")
            utils.plot_results(t1, t2, t, st, mccm, velocity, azimuth, lts_dict, config, array_params, skip_chans)
        except:
            import traceback
            print("Something went wrong making the plot:")
            print(traceback.format_exc())

    #### Write output files ####
    if 'OUT_VALVE_DIR' in config.keys():
        try:
            print('Writing CSV file...')
            t = np.array([utc(dates.num2date(ti)).strftime('%Y-%m-%d %H:%M:%S') for ti in t])
            sta_name = st[0].stats.station
            utils.write_valve_file(t2, t, pressure, azimuth, velocity, mccm, sigma_tau, sta_name, config)
        except:
            import traceback
            print('Something went wrong writing the CSV file:')
            print(traceback.format_exc())

    if 'OUT_ASCII_DIR' in config.keys():
        try:
            print('Writing ASCII file...')
            t = np.array([utc(dates.num2date(ti)).strftime('%Y-%m-%d %H:%M:%S') for ti in t])
            utils.write_ascii_file(t2, t, pressure, azimuth, velocity, mccm, sigma_tau, array_name, config)
        except:
            import traceback
            print('Something went wrong writing the ASCII file:')
            print(traceback.format_exc())

    return


if __name__ == "__main__":
    """
    Main entry point for the script. Handles argument parsing, configuration loading,
    and processing of arrays.
    """
    timer_0 = time.time()

    args = parse_args()  # Parse command-line arguments
    config_file = args.config
    config = utils.load_config(config_file)  # Load configuration

    config["plot"] = False if args.no_plot else True

    T0, delay = get_starttime(config, args)  # Determine start time and delay

    if os.getenv("FROMCRON") == "yep":
        # Set up logging if running from a cron job
        utils.write_to_log(T0.strftime("%Y-%m-%d"), config)

    time.sleep(delay)  # Pause to allow for data latency to catch up
    timer_0 += delay

    print(f"Start time: {utc.utcnow()}")

    if args.array:
        # Process a single array if specified
        timer_tmp = time.time()
        process_array(config, args.array.replace("_", " "), T0)
        dt = time.time() - timer_tmp
        print(f"{dt:.1f} seconds to process {args.array.replace('_', ' ')}")
    else:
        # Process all arrays in the configuration
        for array_name in config["array_list"]:
            timer_tmp = time.time()
            process_array(config, array_name, T0)
            dt = time.time() - timer_tmp
            print(f"{dt:.1f} seconds to process {array_name}\n\n")

    # Write out the new HTML file
    write_html(config)

    print(f"{time.time()-timer_0:.1f} seconds to process all")
    print(f"Finish time: {utc.utcnow()}")
