import os
import sys
import logging
from pathlib import Path
import yaml
import numpy as np
import pandas as pd
import jinja2
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

    default_config_file = Path(__file__).parent.parent.parent / "config/ipensive_config.yml"

    if "IPENSIVE_CONFIG" in os.environ:
        env_config_file = Path(os.environ["IPENSIVE_CONFIG"])
        if env_config_file.exists():
            default_config_file = env_config_file
            my_log.info(f"Using config file from IPENSIVE_CONFIG: {default_config_file}")
        else:
            my_log.warning(f"IPENSIVE_CONFIG does not exist: {env_config_file}")
            raise Exception(f"{env_config_file} does not exist")
    else:
        my_log.info(f"Using default config file: {default_config_file}")

    return default_config_file


def load_array_config(array_file, ipensive_config):
    """
    Load array configuration from a YAML file.

    Args:
        array_file (str): Path to the array configuration file.

    Returns:
        dict: Array configuration dictionary.
    """
    with open(array_file, "r") as file:
        arrays_config = yaml.safe_load(file)

    ###### Load data source configuration ######
    client = get_obspy_client(ipensive_config)

    network_list = list(arrays_config["NETWORKS"].keys())
    array_list = []
    for net in network_list:
        for array in arrays_config["NETWORKS"][net]:
            arrays_config[array]["NETWORK_NAME"] = net
            arrays_config[array]["ARRAY_NAME"] = array
            array_list.append(array)

            if "CLIENT_TYPE" not in arrays_config[array]:
                arrays_config[array]["CLIENT"] = client
            else:
                arrays_config[array]["CLIENT"] = get_obspy_client(arrays_config[array])

            if "STATION_XML" in arrays_config:
                arrays_config[array]["STATION_XML"] = Path(arrays_config["STATION_XML"]).resolve()
            else:
                arrays_config[array]["STATION_XML"] = Path(ipensive_config["STATION_XML"]).resolve()
            if "TARGETS_FILE" in arrays_config:
                arrays_config[array]["TARGETS_FILE"] = Path(arrays_config["TARGETS_FILE"]).resolve()
            else:
                arrays_config[array]["TARGETS_FILE"] = Path(ipensive_config["TARGETS_FILE"]).resolve()

    arrays_config.pop("PARAMS")
    arrays_config.pop("NETWORKS")
    if "STATION_XML" in arrays_config:
        arrays_config.pop("STATION_XML")
    if "TARGETS_FILE" in arrays_config:
        arrays_config.pop("TARGETS_FILE")

    return arrays_config, network_list, array_list


def get_network_dict(config):
    """Get a dictionary mapping network names to their corresponding arrays.

    Args:
        config (dict): Configuration dictionary.

    Returns:
        dict: Dictionary mapping network names to lists of arrays.
    """
    
    network_dict = {}
    for array in config["array_list"]:
        network_name = config[array]["NETWORK_NAME"]
        if network_name not in network_dict:
            network_dict[network_name] = []
        network_dict[network_name].append(array)

    return network_dict


def load_ipensive_config(config_file):
    """
    Load configuration from a YAML file and initialize network and array settings.

    Args:
        config_file (str): Path to the configuration file.

    Returns:
        dict: Configuration dictionary with additional metadata.
        pathlib.Path: Path to the arrays configuration file.
    """

    ###### Load main iPensive config ######
    with open(config_file, "r") as file:
        ipensive_config = yaml.safe_load(file)

    ###### Load array configurations ######
    if "ARRAYS_CONFIG" in os.environ:
        array_file = [os.environ["ARRAYS_CONFIG"]]
    elif "ARRAYS_CONFIG" in ipensive_config:
        array_file = ipensive_config["ARRAYS_CONFIG"]
    else:
        array_file = Path(__file__).parent.parent.parent / "config/arrays_config.yml"
    if not isinstance(array_file, list):
        array_file = [array_file]
    array_file = [Path(af).resolve() for af in array_file]

    for i, f in enumerate(array_file):
        my_log.info(f"Loading arrays config file: {f}")

        if i == 0:
            ARRAYS_CONFIG, net_list, array_list = load_array_config(f, ipensive_config)
        if i > 0:
            tmp_arrays_config, tmp_net_list, tmp_array_list = load_array_config(f, ipensive_config)
            ARRAYS_CONFIG.update(tmp_arrays_config)
            net_list.extend(tmp_net_list)
            array_list.extend(tmp_array_list)

    ##### Add network and array summary information
    ARRAYS_CONFIG["array_list"] = array_list

    if "EXTRA_LINKS" in ipensive_config:
        ARRAYS_CONFIG["EXTRA_LINKS"] = ipensive_config["EXTRA_LINKS"]
    else:
        ARRAYS_CONFIG["EXTRA_LINKS"] = []

    if "LATENCY" in ipensive_config:
        ARRAYS_CONFIG["LATENCY"] = ipensive_config["LATENCY"]
    else:
        ARRAYS_CONFIG["LATENCY"] = 0  # default latency in seconds

    ###### Load data output configuration ######
    ARRAYS_CONFIG["OUT_WEB_DIR"] = Path(ipensive_config["OUT_WEB_DIR"])
    ARRAYS_CONFIG["LOGS_DIR"] = Path(ipensive_config["LOGS_DIR"])
    if "OUT_ASCII_DIR" in ipensive_config and ipensive_config["OUT_ASCII_DIR"]:
        ARRAYS_CONFIG["OUT_ASCII_DIR"] = Path(ipensive_config["OUT_ASCII_DIR"])
    if "OUT_VALVE_DIR" in ipensive_config and ipensive_config["OUT_VALVE_DIR"]:
        ARRAYS_CONFIG["OUT_VALVE_DIR"] = Path(ipensive_config["OUT_VALVE_DIR"])

    ARRAYS_CONFIG["array_config_files"] = array_file

    ARRAYS_CONFIG["DURATION"] = 600  # default duration in seconds
    for array in array_list:
        ARRAYS_CONFIG[array]["DURATION"] = ARRAYS_CONFIG["DURATION"]  # default duration in seconds

    ARRAYS_CONFIG["NETWORKS"] = get_network_dict(ARRAYS_CONFIG)

    return ARRAYS_CONFIG


def get_pngfile_path(t, array_name, config):
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


def get_target_backazimuth(st, array_params):
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
    DF = pd.read_csv(array_params["TARGETS_FILE"])  # Load target locations from a CSV file

    tmp_targets = []
    tmp_baz = []
    
    if isinstance(array_params["TARGETS"], dict):
        for target in array_params["TARGETS"]:
            tmp_targets.append(target)
            tmp_baz.append(array_params["TARGETS"][target])
    elif isinstance(array_params["TARGETS"], list):
        for target in array_params["TARGETS"]:
            # If the target is a string, look up its coordinates in the CSV
            df = DF[DF["Target"] == target]
            _, baz, _ = gps2dist_azimuth(lat0, lon0, df.iloc[0]["Latitude"], df.iloc[0]["Longitude"])
            tmp_targets.append(target)
            tmp_baz.append(baz)

    # Add the calculated backazimuths to the array parameters
    for t, baz in zip(tmp_targets, tmp_baz):
        array_params[t] = baz

    return array_params


def write_html(config):
    """
    Generate an HTML file for web output.

    Args:
        config (dict): Configuration dictionary.

    Returns:
        None
    """
    
    template_file = Path(__file__).parent.parent.parent / "templates" / "index.template"
    with open(template_file, "r") as f:
        template = jinja2.Template(f.read())
    html = template.render(
        networks=list(config['NETWORKS'].keys()),
        arrays=config["NETWORKS"],
        extra_links=config["EXTRA_LINKS"],
    )
    out_file = config["OUT_WEB_DIR"] / "index.html"
    with open(out_file, "w") as f:
        f.write(html)


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


def write_data_files(t2, st, df, config):
    """Write data files for the given time, stream, and dataframe.

    Args:
        t2 (obspy.UTCDateTime): Current time.
        st (obspy.Stream): Stream containing seismic traces.
        df (pd.DataFrame): DataFrame containing metadata and results.
        config (dict): Configuration dictionary.
    """
    
    df["Time"] = [dates.num2date(ti).strftime('%Y-%m-%d %H:%M:%S') for ti in df.Time]
    if "OUT_VALVE_DIR" in config and config["OUT_VALVE_DIR"]:
        try:
            sta_name = st[0].stats.station
            write_valve_file(t2, df, sta_name, config)
        except Exception:
            import traceback
            my_log.error("Something went wrong writing the Valve CSV file:")
            my_log.error(traceback.format_exc())

    if "OUT_ASCII_DIR" in config and config["OUT_ASCII_DIR"]:
        try:
            write_ascii_file(t2, df, config)
        except Exception:
            import traceback
            my_log.error("Something went wrong writing the ASCII data file:")
            my_log.error(traceback.format_exc())


def write_ascii_file(t2, tmp_df, config):
    """
    Write results to an ASCII file.

    Args:
        t2 (obspy.UTCDateTime): Current time.
        df (pd.Dataframe): Dataframe containing:
                            Time: List of timestamps.
                            Array: Array name.
                            Pressure: Pressure values.
                            Azimuth: Azimuth values.
                            Velocity: Velocity values.
                            MCCM: MCCM values.
                            Sigma_tau: Sigma_tau values.
                            Vel_err: Confidence interval for velocity.
                            Baz_err: Confidence interval for backazimuth.
        name (str): Array name.
        config (dict): Configuration dictionary.

    Returns:
        None
    """

    my_log.info("Writing ASCII file...")

    array_name = tmp_df.iloc[0]["Array"]
    t1 = t2 - config[array_name]["DURATION"]
    array_path_name = array_name.replace(" ", "_")

    # Create output directories
    array_dir = config["OUT_ASCII_DIR"] / array_path_name
    month_dir = array_dir / t1.strftime("%Y-%m")
    month_dir.mkdir(parents=True, exist_ok=True)

    filename = month_dir / f"{array_path_name}_{t1.strftime('%Y-%m-%d')}.txt"

    # Adjust azimuth values to be within 0-360 degrees
    tmp_df.loc[tmp_df["Azimuth"] < 0, "Azimuth"] += 360

    # Convert time to datetime
    tmp_df["Time"] = pd.to_datetime(tmp_df["Time"])
    tmp_df = tmp_df[tmp_df["Time"] <= t2.strftime('%Y-%m-%d %H:%M:%S')]

    # Append to or overwrite the existing file
    if filename.exists():
        df = pd.read_csv(filename, sep="\t", parse_dates=["Time"])
        df = df[(df["Time"] <= t1.strftime("%Y-%m-%d %H:%M:%S")) | (df["Time"] > t2.strftime("%Y-%m-%d %H:%M:%S"))]
        df = df.rename(columns={"rms": "Sigma_tau"})
        df = pd.concat([df, tmp_df])
        df = df.sort_values("Time")
    else:
        df = tmp_df

    # Round values for better readability
    df = df.round({"Azimuth": 1, "Velocity": 1, "MCCM": 2, "Pressure": 3, "Sigma_tau": 2, "Vel_err": 1, "Baz_err": 1})

    # Save the DataFrame to a file
    df.to_csv(filename, index=False, header=True, sep="\t")

    return


def write_valve_file(t2, df, name, config):
    """
    Write results to a CSV file for valve output.

    Args:
        t2 (obspy.UTCDateTime): Current time.
        df (pd.Dataframe): Dataframe containing:
                            Time: List of timestamps.
                            Array: Array name.
                            Pressure: Pressure values.
                            Azimuth: Azimuth values.
                            Velocity: Velocity values.
                            MCCM: MCCM values.
                            Sigma_tau: Sigma_tau values.
                            Vel_err: Confidence interval for velocity.
                            Baz_err: Confidence interval for backazimuth.                            
        config (dict): Configuration dictionary.

    Returns:
        None
    """

    my_log.info('Writing CSV file for Valve...')

    # Create a DataFrame with the results
    A = df.rename(columns={"Time": "TIMESTAMP", "Array": "CHANNEL"})

    A = A[["TIMESTAMP", "CHANNEL", "Azimuth", "Velocity", "MCCM", "Pressure", "Sigma_tau", "Vel_err", "Baz_err"]]

    # Save the DataFrame to a CSV file
    out_valve_dir = config["OUT_VALVE_DIR"]
    out_valve_dir.mkdir(parents=True, exist_ok=True)
    filename = out_valve_dir / f"{name}_{t2.strftime('%Y%m%d-%H%M')}.txt"
    A.to_csv(filename, index=False, header=True, float_format="%.3f")

    return
