#!/home/rtem/miniconda3/envs/ipensive/bin/python
# -*- coding: utf-8 -*-
"""
Created on 24-Apr-2018
Modified on 30-Mar-2022
@author: awech
"""

import os
import sys
sys.dont_write_bytecode = True  # don't write .pyc files (probably slightly faster without this, but more cluttered)
import numpy as np
import jinja2
from obspy.core import UTCDateTime as utc
from matplotlib import dates
import time
import ipensive_utils as utils
import argparse


def parse_args():
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
        help="utc time stamp:YYYYMMDDHHMM (optional, otherwise grabs current utc time)",
    )
    parser.add_argument(
        "-a",
        "--array",
        type=str,
        help="Array name if you want process single array. (Use _ instead of spaces, if necessary)",
    )

    return parser.parse_args()


def get_starttime(config, args):

    date_fmt = "%Y-%m-%d %H:%M"

    if args.time is None:
        T0 = utc.utcnow()  # no time given, use current timestamp
        delay = config["PARAMS"]["LATENCY"] + config["PARAMS"]["WINDOW_LENGTH"]
    else:
        T0 = utc(args.time)
        delay = 0

    T0 = utc(T0.strftime(date_fmt)[:-1] + "0")  # round down to the nearest 10-minute

    return T0, delay


def write_html(config):
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

    params = config[array_name]
    t1 = T0 - params["DURATION"]
    t2 = T0
    fmt = "%Y-%b-%d %H:%M"
    print(f"Processing {array_name}: {t1.strftime(fmt)} - {t2.strftime(fmt)}")

    T1 = t1 - params["TAPER"]
    T2 = t2 + params["WINDOW_LENGTH"] + params["TAPER"]

    print("--- " + params["ARRAY_NAME"] + " ---")
    if os.getenv("FROMCRON") == "yep":
        time.sleep(params["EXTRA_PAUSE"])

    #### download data ####
    if len(params["SCNL"]) < params["MIN_CHAN"]:
        print("Not enough channels defined.")
        return

    st = utils.grab_data(
        params["SCNL"],
        T1,
        T2,
        hostname=params["HOSTNAME"],
        port=params["PORT"],
        fill_value=0,
    )
    # st = utils.add_coordinate_info(st, scnl)
    st = utils.add_metadata(st, config["STATION_XML"])
    params = utils.get_volcano_backazimuth(st, config, params)
    ########################

    #### check for enough data ####
    for tr in st:
        if np.sum(np.abs(tr.data)) == 0:
            st.remove(tr)
    if len(st) < params["MIN_CHAN"]:
        print("Too many blank traces. Skipping.")
        return
    ########################

    #### check for gappy data ####
    for tr in st:
        if np.any([np.any(tr.data == 0)]):
            st.remove(tr)
    if len(st) < params["MIN_CHAN"]:
        print("Too gappy. Skipping.")
        return
    ########################

    #### preprocess data ####
    for tr in st:
        tr.remove_sensitivity(tr.inventory)
    st.detrend("demean")
    for tr in st:
        if tr.stats["sampling_rate"] == 100.0:
            tr.decimate(2)
        if tr.stats["sampling_rate"] != 50.0:
            tr.resample(50.0)
        if tr.stats["sampling_rate"] == 50.0:
            if float(params["FREQMAX"]) < 50 / 4.0:
                tr.decimate(2)
    st.taper(max_percentage=None, max_length=params["TAPER"])
    st.filter(
        "bandpass",
        freqmin=params["FREQMIN"],
        freqmax=params["FREQMAX"],
        corners=2,
        zerophase=True,
    )
    st.trim(t1, t2 + params["WINDOW_LENGTH"])
    ########################

    velocity = []
    azimuth = []
    mccm = []
    t = []
    rms = []
    pressure = []
    for st_win in st.slide(
        window_length=params["WINDOW_LENGTH"],
        step=params["WINDOW_LENGTH"] - params["OVERLAP"],
    ):
        try:
            
            vel, az, rms0, cmax, pk_press = utils.inversion(st_win)
            velocity.append(vel)
            azimuth.append(az)
            mccm.append(np.median(cmax))
            t.append(
                dates.date2num(
                    (st_win[0].stats.starttime + params["WINDOW_LENGTH"] / 2.0).datetime
                )
            )
            rms.append(rms0)
            pressure.append(pk_press)
        except:
            print("Something went wrong in the inversion...")
            continue
    t = np.array(t)
    mccm = np.array(mccm)
    velocity = np.array(velocity)
    azimuth = np.array(azimuth)
    rms = np.array(rms)
    pressure = np.array(pressure)


    try:
        print("Setting up web output folders")
        utils.web_folders(t2, config, params)
        print("Making plot...")
        utils.plot_results(t1, t2, t, st, mccm, velocity, azimuth, config, params)
    except:
        import traceback

        b = traceback.format_exc()
        message = "".join(f"{a}\n" for a in b.splitlines())
        print("Something went wrong making the plot:")
        print(message)

    if 'OUT_VALVE_DIR' in dir(config):
        try:
            print('Writing csv file...')
            t=np.array([utc(dates.num2date(ti)).strftime('%Y-%m-%d %H:%M:%S') for ti in t])
            name=st[0].stats.station
            utils.write_valve_file(t2, t, pressure, azimuth, velocity, mccm, rms, name)
        except:
            import traceback
            b=traceback.format_exc()
            message = ''.join('{}\n'.format(a) for a in b.splitlines())
            print('Something went wrong writing the csv file:')
            print(message)

    if 'OUT_ASCII_DIR' in dir(config):
        try:
            print('Writing csv file...')    
            t=np.array([utc(dates.num2date(ti)).strftime('%Y-%m-%d %H:%M:%S') for ti in t])
            name=array['Name']
            utils.write_ascii_file(t2, t, pressure, azimuth, velocity, mccm, rms, name)
        except:
            import traceback
            b=traceback.format_exc()
            message = ''.join('{}\n'.format(a) for a in b.splitlines())
            print('Something went wrong writing the csv file:')
            print(message)

    return


if __name__ == "__main__":

    timer_0 = time.time()

    args = parse_args()
    config_file = args.config
    config = utils.load_config(config_file)
    T0, delay = get_starttime(config, args)

    if os.getenv("FROMCRON") == "yep":
        # Set up logging
        utils.write_to_log(T0.strftime("%Y-%m-%d"), config)

    time.sleep(delay)  # Pause to allow for data latency to catch up
    timer_0 += delay

    print(f"Start time: {utc.utcnow()}")

    if args.array:
        timer_tmp = time.time()
        process_array(config, args.array.replace("_", " "), T0)
        dt = time.time() - timer_tmp
        print(f"{dt:.1f} seconds to process {args.array.replace("_", " ")}")
    else:
        for array_name in config["array_list"]:
            timer_tmp = time.time()
            process_array(config, array_name, T0)
            dt = time.time() - timer_tmp
            print(f"{dt:.1f} seconds to process {array_name}")

    # Write out the new HTML file
    write_html(config)

    print(f"{time.time()-timer_0:.1f} seconds to process all")
    print(f"Finish time: {utc.utcnow()}")
