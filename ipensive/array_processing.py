"""
Created on 24-Apr-2018
Modified on 15-May-2025
@author: awech
"""

import os
import argparse
import logging
import time
from pathlib import Path
import numpy as np
import jinja2
from obspy import UTCDateTime as utc
from . import ipensive_utils as utils
from .plotting_utils import plot_results
from .data_utils import grab_data, QC_data, preprocess_data
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
        delay = 0
    else:
        # Use the current UTC time and add latency and window length
        T0 = utc.utcnow()
        delay = config["PARAMS"]["LATENCY"] + config["PARAMS"]["WINDOW_LENGTH"]
        my_log.info(f"Waiting {delay:g} seconds")

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
    
    template_file = Path(__file__).parent.parent / "templates" / "index.template"
    with open(template_file, "r") as f:
        template = jinja2.Template(f.read())
    html = template.render(
        networks=config["network_list"],
        arrays=config["NETWORKS"],
        extra_links=config["EXTRA_LINKS"],
    )
    out_file = Path(config["OUT_WEB_DIR"]) / "index.html"
    with open(out_file, "w") as f:
        f.write(html)


def process_array(config, array_name, T0):
    """
    Process a single array for the specified time window.

    Args:
        config (dict): Configuration dictionary.
        array_name (str): Name of the array to process.
        T0 (obspy.UTCDateTime): End time of the processing window.

    Returns:
        None
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
        return
    # Check for latency delay
    if os.getenv("FROMCRON") == "yep":
        if "EXTRA_PAUSE" in array_params:
            time.sleep(array_params["EXTRA_PAUSE"])  # Pause if running from a cron job

    # Download data
    st = grab_data(
        array_params["CLIENT"],
        array_params["NSLC"],
        T1,
        T2,
    )

    # Check data quality
    good_data, skip_chans = QC_data(st, array_params)
    if not good_data:
        return

    # Add metadata or coordinates
    st = add_metadata(st, config, array_name, skip_chans)

    # Add volcano backazimuths
    array_params = utils.get_target_backazimuth(st, config, array_params)

    # Preprocess data
    st = preprocess_data(st, t1, t2, skip_chans, array_params)


    # Perform array processing
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


    # Generate plots
    if config["plot"]:
        try:
            utils.web_folders(t2, config, array_params)
            for plotsize in ["big", "small"]:
                plot_results(t1, t2, t, st, mccm, velocity, azimuth, lts_dict, skip_chans, config, array_params, plotsize)
        except Exception:
            import traceback
            my_log.error("Something went wrong making the plot:")
            my_log.error(traceback.format_exc())


    # Write output files
    if 'OUT_VALVE_DIR' in config.keys():
        try:
            sta_name = st[0].stats.station
            utils.write_valve_file(t2, t, pressure, azimuth, velocity, mccm, sigma_tau, sta_name, config)
        except Exception:
            import traceback
            my_log.error('Something went wrong writing the CSV file:')
            my_log.error(traceback.format_exc())

    if 'OUT_ASCII_DIR' in config.keys() and config['OUT_ASCII_DIR']:
        try:
            utils.write_ascii_file(t2, t, pressure, azimuth, velocity, mccm, sigma_tau, array_name, config)
        except Exception:
            import traceback
            my_log.error('Something went wrong writing the ASCII file:')
            my_log.error(traceback.format_exc())

    return
