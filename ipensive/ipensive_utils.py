import os
import sys
import logging
from pathlib import Path
import yaml
import numpy as np
import pandas as pd
from matplotlib import dates
from obspy import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth
from .data_utils import get_obspy_client

my_log = logging.getLogger(__name__)


class StreamToLogger(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    """
    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())

    def flush(self):
        pass


def get_config_file():
    """
    Get the path to the configuration file.

    Returns:
        pathlib.Path: Path to the configuration file.
    """

    default_config_file = Path(__file__).parent.parent / "config" / "ipensive_config.yml"

    if "IPENSIVE_CONFIG" in os.environ:
        env_config_file = Path(os.environ["IPENSIVE_CONFIG"])
        if env_config_file.exists():
            default_config_file = env_config_file
            print(f"Using config file from IPENSIVE_CONFIG: {default_config_file}")
            my_log.info(f"Using config file from IPENSIVE_CONFIG: {default_config_file}")
        else:
            print(f"IPENSIVE_CONFIG does not exist: {env_config_file}")
            my_log.warning(f"IPENSIVE_CONFIG does not exist: {env_config_file}")
            raise Exception(f"{env_config_file} does not exist")
    else:
        print(f"Using default config file: {default_config_file}")
        my_log.info(f"Using default config file: {default_config_file}")

    return default_config_file


def load_config(config_file):
    """
    Load configuration from a YAML file and initialize network and array settings.

    Args:
        config_file (str): Path to the configuration file.

    Returns:
        dict: Configuration dictionary with additional metadata.
    """

    ###### Load main iPensive config ######
    with open(config_file, "r") as file:
        ipensive_config = yaml.safe_load(file)


    ###### Load array configurations ######
    if "ARRAYS_CONFIG" in os.environ:
        array_file = os.environ["ARRAYS_CONFIG"]
    elif "ARRAYS_CONFIG" in ipensive_config.keys():
        array_file = ipensive_config["ARRAYS_CONFIG"]
    else:
        array_file = Path(__file__).parent.parent / "config" / "arrays_config.yml"
    print(f"Using arrays config file: {array_file}")
    my_log.info(f"Using arrays config file: {array_file}")
    with open(array_file, "r") as file:
        arrays_config = yaml.safe_load(file)


    ###### Load target & metadata info ######
    arrays_config["STATION_XML"] = Path(ipensive_config["STATION_XML"])
    arrays_config["TARGETS_FILE"] = Path(ipensive_config["TARGETS_FILE"])
    if "EXTRA_LINKS" not in arrays_config.keys():
        arrays_config["EXTRA_LINKS"] = []


    ##### Extract network and array information
    all_nets = list(arrays_config["NETWORKS"].keys())
    array_list = []
    for net in all_nets:
        for array in arrays_config["NETWORKS"][net]:
            arrays_config[array]["NETWORK_NAME"] = net
            arrays_config[array]["ARRAY_NAME"] = array
            array_list.append(array)
    arrays_config["network_list"] = all_nets
    arrays_config["array_list"] = array_list


    ###### Load data output configuration ######
    arrays_config["OUT_WEB_DIR"] = Path(ipensive_config["OUT_WEB_DIR"])
    arrays_config["OUT_ASCII_DIR"] = Path(ipensive_config["OUT_ASCII_DIR"])
    arrays_config["LOGS_DIR"] = Path(ipensive_config["LOGS_DIR"])

    ###### Load data source configuration ######
    client = get_obspy_client(ipensive_config)
    
    # Assign or update the client for each array
    for array in arrays_config["array_list"]:
        if "CLIENT_TYPE" not in arrays_config[array]:
            arrays_config[array]["CLIENT"] = client
        else:
            arrays_config[array]["CLIENT"] = get_obspy_client(arrays_config[array])

    return arrays_config


def get_file_path(t, array_name, config):
    """
    Get the file path for .png output of a specific array and time.

    Args:
        t (pandas Timestamp): The time for which to get the file path.
        array_name (str): The name of the array.
        config (dict): The configuration dictionary.

    Returns:
        pathlib.Path: The file path for the .png file of the specified array and time.
    """

    array_dict = config.get(array_name, {})
    network_dir = config["OUT_WEB_DIR"] / array_dict["NETWORK_NAME"]
    array_dir = network_dir / array_dict["ARRAY_NAME"]
    year_dir = array_dir / str(t.year)
    julian_day_dir = year_dir / str(t.day_of_year)

    time = t.strftime("%Y%m%d-%H%M")
    file = julian_day_dir / f"{array_name}_{time}.png"

    return file


def setup_logging(day, config, arg_opt=None):
    """
    Write logs to a file for a specific day.

    Args:
        day (str): Date string in UTC.
        config (dict): Configuration dictionary.
    """

    # Move path to `load_config`
    # Determine the logs directory
    if 'LOGS_DIR' in config:
        if config["LOGS_DIR"] is not None:
            logs_dir = Path(config["LOGS_DIR"])
        else:
            logs_dir = None
    else:
        logs_dir = Path(__file__).parent / 'logs'

    if logs_dir is not None and os.getenv("FROMCRON") == "yep":
        # Create year and month directories if they don't exist
        year = UTCDateTime(day).strftime('%Y')
        month = UTCDateTime(day).strftime('%Y-%m')
        year_dir = logs_dir / year
        month_dir = year_dir / month

        year_dir.mkdir(parents=True, exist_ok=True)
        month_dir.mkdir(parents=True, exist_ok=True)

        # Create the log file
        log_file = month_dir / f"{UTCDateTime(day).strftime('%Y-%m-%d')}.log"
    elif arg_opt:
        log_file = arg_opt
    else:
        log_file = None


    if log_file:
        logging.basicConfig(
            filename=log_file,
            filemode="a",
            level=logging.INFO,
            format="%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
    else:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
    my_log = logging.getLogger(__name__)
    sys.stdout = StreamToLogger(my_log, logging.INFO)
    sys.stderr = StreamToLogger(my_log, logging.ERROR)


def get_target_backazimuth(st, config, array_params):
    """
    Calculate the backazimuths for target locations relative to the array's average coordinates.

    Args:
        st (Stream): ObsPy Stream object containing traces with coordinate metadata.
        config (dict): Configuration dictionary.
        array_params (dict): Parameters for the array, including target information.

    Returns:
        dict: Updated array parameters with backazimuths for each target.
    """
    # Calculate the average latitude and longitude of the array
    lon0 = np.nanmean([tr.stats.coordinates.longitude for tr in st])
    lat0 = np.nanmean([tr.stats.coordinates.latitude for tr in st])
    DF = pd.read_csv(config["TARGETS_FILE"])  # Load target locations from a CSV file

    tmp_targets = []
    tmp_baz = []
    for target in array_params["TARGETS"]:
        if type(target) is str:
            # If the target is a string, look up its coordinates in the CSV
            df = DF[DF["Target"] == target]
            _, baz, _ = gps2dist_azimuth(lat0, lon0, df.iloc[0]["Latitude"], df.iloc[0]["Longitude"])
            tmp_targets.append(target)
            tmp_baz.append(baz)
        elif type(target) is dict:
            # If the target is a dictionary, use the provided backazimuth
            t_name = list(target.keys())[0]
            tmp_targets.append(t_name)
            tmp_baz.append(target[t_name])

    # Add the calculated backazimuths to the array parameters
    for t, baz in zip(tmp_targets, tmp_baz):
        array_params[t] = baz

    return array_params


def web_folders(t2, config, params):
    """
    Create directory structure for web output.

    Args:
        t2 (obspy.UTCDateTime): Current time.
        config (dict): Configuration dictionary.
        params (dict): Array parameters.

    Returns:
        None
    """

    my_log.info("Setting up web output folders")

    # Create the base output directory if it doesn't exist
    out_web_dir = config["OUT_WEB_DIR"]
    out_web_dir.mkdir(parents=True, exist_ok=True)

    # Create subdirectories for the network, array, year, and Julian day
    network_dir = out_web_dir / params["NETWORK_NAME"]
    array_dir = network_dir / params["ARRAY_NAME"]
    year_dir = array_dir / str(t2.year)
    julian_day_dir = year_dir / f"{t2.julday:03d}"

    julian_day_dir.mkdir(parents=True, exist_ok=True)
    return


def write_ascii_file(t2, t, pressure, azimuth, velocity, mccm, rms, name, config):
    """
    Write results to an ASCII file.

    Args:
        t2 (obspy.UTCDateTime): Current time.
        t (list): List of timestamps.
        pressure (list): Pressure values.
        azimuth (list): Azimuth values.
        velocity (list): Velocity values.
        mccm (list): MCCM values.
        rms (list): RMS values.
        name (str): Array name.
        config (dict): Configuration dictionary.

    Returns:
        None
    """

    my_log.info('Writing ASCII file...')
    
    t1 = t2 - config[name]["DURATION"]
    tmp_name = name.replace(' ', '_')

    # Create output directories
    array_dir = config["OUT_ASCII_DIR"] / tmp_name
    month_dir = array_dir / t1.strftime('%Y-%m')
    month_dir.mkdir(parents=True, exist_ok=True)

    filename = month_dir / f"{tmp_name}_{t1.strftime('%Y-%m-%d')}.txt"

    # Adjust azimuth values to be within 0-360 degrees
    azimuth[azimuth < 0] += 360

    t = np.array([UTCDateTime(dates.num2date(ti)).strftime('%Y-%m-%d %H:%M:%S') for ti in t])

    # Create a DataFrame with the results
    tmp = pd.DataFrame({
        'Time': t,
        'Array': name,
        'Azimuth': azimuth,
        'Velocity': 1000 * velocity,  # Convert velocity to m/s
        'MCCM': mccm,
        'Pressure': pressure,
        'rms': rms
    })
    tmp['Time'] = pd.to_datetime(tmp['Time'])
    tmp = tmp[tmp['Time'] <= t2.strftime('%Y-%m-%d %H:%M:%S')]

    # Append to or overwrite the existing file
    if filename.exists():
        df = pd.read_csv(filename, sep='\t', parse_dates=['Time'])
        df = df[(df['Time'] <= t1.strftime('%Y-%m-%d %H:%M:%S')) | (df['Time'] > t2.strftime('%Y-%m-%d %H:%M:%S'))]
        df = pd.concat([df, tmp])
        df = df.sort_values('Time')
    else:
        df = tmp

    # Round values for better readability
    df = df.round({'Azimuth': 1, 'Velocity': 1, 'MCCM': 2, 'Pressure': 3, 'rms': 1})

    # Save the DataFrame to a file
    df.to_csv(filename, index=False, header=True, sep='\t')
    return


def write_valve_file(t2, t, pressure, azimuth, velocity, mccm, rms, name, config):
    """
    Write results to a CSV file for valve output.

    Args:
        t2 (obspy.UTCDateTime): Current time.
        t (list): List of timestamps.
        pressure (list): Pressure values.
        azimuth (list): Azimuth values.
        velocity (list): Velocity values.
        mccm (list): MCCM values.
        rms (list): RMS values.
        name (str): Array name.
        config (dict): Configuration dictionary.

    Returns:
        None
    """

    my_log.info('Writing CSV file for Valve...')

    t = np.array([UTCDateTime(dates.num2date(ti)).strftime('%Y-%m-%d %H:%M:%S') for ti in t])

    # Create a DataFrame with the results
    A = pd.DataFrame({
        'TIMESTAMP': t,
        'CHANNEL': name,
        'Azimuth': azimuth,
        'Velocity': 1000 * velocity,  # Convert velocity to m/s
        'MCCM': mccm,
        'Pressure': pressure,
        'rms': rms
    })

    # Save the DataFrame to a CSV file
    out_valve_dir = config["OUT_VALVE_DIR"]
    out_valve_dir.mkdir(parents=True, exist_ok=True)
    filename = out_valve_dir / f"{name}_{t2.strftime('%Y%m%d-%H%M')}.txt"
    A.to_csv(filename, index=False, header=True, sep=',', float_format='%.3f')
    return

