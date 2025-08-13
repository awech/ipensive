import os
import sys
import logging
from pathlib import Path
from time import sleep
import yaml
import numpy as np
import pandas as pd
from obspy import Stream, UTCDateTime, read_inventory
from obspy.clients.earthworm import Client as EWClient
from obspy.clients.fdsn import Client as FDSNClient
from obspy.clients.filesystem.sds import Client as SDSClient
from obspy.clients.seedlink import Client as SLClient
from obspy.geodetics.base import gps2dist_azimuth
from obspy.core.util import AttribDict

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
        Path: Path to the configuration file.
    """

    default_config = Path(__file__).parent.parent / "config" / "ipensive_config.yml"

    if "IPENSIVE_CONFIG" in os.environ:
        env_config_file = Path(os.environ["IPENSIVE_CONFIG"])
        if env_config_file.exists():
            default_config = env_config_file
            print(f"Using config file from IPENSIVE_CONFIG: {default_config}")
            my_log.info(f"Using config file from IPENSIVE_CONFIG: {default_config}")
        else:
            print(f"IPENSIVE_CONFIG does not exist: {env_config_file}")
            my_log.warning(f"IPENSIVE_CONFIG does not exist: {env_config_file}")
            raise Exception(f"{env_config_file} does not exist")
    else:
        print(f"Using default config file: {default_config}")
        my_log.info(f"Using default config file: {default_config}")

    return default_config


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
        config = yaml.safe_load(file)


    ###### Load target & metadata info ######
    config["STATION_XML"] = ipensive_config["STATION_XML"]
    config["TARGETS_FILE"] = ipensive_config["TARGETS_FILE"]
    if "EXTRA_LINKS" not in config.keys():
        config["EXTRA_LINKS"] = []


    ##### Extract network and array information
    all_nets = list(config["NETWORKS"].keys())
    array_list = []
    for net in all_nets:
        for array in config["NETWORKS"][net]:
            config[array]["NETWORK_NAME"] = net
            config[array]["ARRAY_NAME"] = array
            array_list.append(array)
    config["network_list"] = all_nets
    config["array_list"] = array_list


    ###### Load data output configuration ######
    config["OUT_WEB_DIR"] = ipensive_config["OUT_WEB_DIR"]
    config["OUT_ASCII_DIR"] = ipensive_config["OUT_ASCII_DIR"]
    config["LOGS_DIR"] = ipensive_config["LOGS_DIR"]

    ###### Load data source configuration ######
    client = get_obspy_client(ipensive_config)
    
    # Assign or update the client for each array
    for array in config["array_list"]:
        if "CLIENT_TYPE" not in config[array].keys():
            config[array]["CLIENT"] = client
        else:
            config[array]["CLIENT"] = get_obspy_client(config[array])

    return config


def get_obspy_client(config):
    """
    Initialize an ObsPy client based on the configuration.

    Args:
        config (dict): Configuration dictionary for the client.

    Returns:
        ObsPy client object.
    """

    if "TIMEOUT" not in config.keys():
        config["TIMEOUT"] = 30

    if config["CLIENT_TYPE"].lower() == "fdsn":
        client = FDSNClient(config["HOSTNAME"], timeout=config["TIMEOUT"])
        client.name = config["HOSTNAME"]

    elif config["CLIENT_TYPE"].lower() == "local_fdsn":
        client = FDSNClient("IRIS", service_mappings={"dataselect": config["LOCAL_FDSN"]}, timeout=config["TIMEOUT"])
        client.name = config["LOCAL_FDSN"]

    elif config["CLIENT_TYPE"].lower() == "sds":
        client = SDSClient(config["DIRECTORY"], timeout=config["TIMEOUT"])
        if "FMT" in list(config.keys()):
            client.FMTSTR = config["FMT"]
        client.name = config["DIRECTORY"]

    elif config["CLIENT_TYPE"].lower() == "earthworm":
        client = EWClient(config["HOSTNAME"], config["PORT"], timeout=config["TIMEOUT"])
        client.name = config["HOSTNAME"]

    elif config["CLIENT_TYPE"].lower() == "seedlink":
        if "PORT" in list(config.keys()):
            client = SLClient(config["HOSTNAME"], config["PORT"], timeout=config["TIMEOUT"])
        else:
            client = SLClient(config["HOSTNAME"], timeout=config["TIMEOUT"])
        client.name = config["HOSTNAME"]
        
    return client


def get_file_path(t, array_name, config):

    out_web_dir = Path(config["OUT_WEB_DIR"])
    array_dict = config.get(array_name, {})
    network_dir = out_web_dir / array_dict["NETWORK_NAME"]
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
    

def check_inventory(tr, inv):
    """
    Check if a trace exists in the inventory.

    Args:
        tr (Trace): ObsPy Trace object.
        inv (Inventory): ObsPy Inventory object.

    Returns:
        bool: True if the trace exists in the inventory, False otherwise.
    """
    inv_test = inv.select(
        network=tr.stats.network,
        station=tr.stats.station,
        location=tr.stats.location,
        channel=tr.stats.channel,
        starttime=tr.stats.starttime,
        endtime=tr.stats.starttime,
    )
    value = True if len(inv_test) > 0 else False
    return value


def check_FDSN(tr, client):
    """
    Check if a trace exists in the FDSN client.

    Args:
        tr (Trace): ObsPy Trace object.
        client (FDSNClient): ObsPy FDSN client.

    Returns:
        bool: True if the trace exists in the FDSN client, False otherwise.
    """
    value = True
    try:
        client.get_stations(
            network=tr.stats.network,
            station=tr.stats.station,
            location=tr.stats.location,
            channel=tr.stats.channel,
            starttime=tr.stats.starttime,
            endtime=tr.stats.starttime,
            level="response",
        )
    except Exception as err:
        if "No data available for request." in err.args[0]:
            value = False
    return value


def add_coordinate_info(st, config, array_name):
    """
    Add coordinate information to traces in a stream.

    Args:
        st (Stream): ObsPy Stream object.
        config (dict): Configuration dictionary.
        array_name (str): Name of the array.

    Returns:
        Stream: Stream with updated coordinate information.
    """
    array_params = config[array_name]
    nslc_params = array_params["NSLC"]

    for tr in st:
        tmp_lat = nslc_params[tr.id.replace("--", "")]["lat"]
        tmp_lon = nslc_params[tr.id.replace("--", "")]["lon"]
        tr.stats.coordinates = AttribDict({
            'latitude': tmp_lat,
            'longitude': tmp_lon,
            'elevation': 0.0
        })
    return st


def add_metadata(st, config, skip_chans=[]):
    """
    Add metadata to traces in a stream.

    Args:
        st (Stream): ObsPy Stream object.
        config (dict): Configuration dictionary.
        skip_chans (list): List of channels to skip.

    Returns:
        Stream: Stream with updated metadata.
    """
    import warnings

    warnings.simplefilter("ignore", UserWarning, append=True)
    if "STATION_XML" in config.keys():
        inventory = read_inventory(config["STATION_XML"])

    for tr in st:
        if tr.id in skip_chans:
            my_log.warning(f"Skipping metadata for {tr.id} due to missing data")
            tr.stats.coordinates = AttribDict({
                'latitude': np.nan,
                'longitude': np.nan,
                'elevation': 0.0
            })
            continue
        my_log.info(f"Getting metadata for {tr.id}")
        if check_inventory(tr, inventory):
            inv = inventory.select(
                network=tr.stats.network,
                station=tr.stats.station,
                location=tr.stats.location,
                channel=tr.stats.channel,
                starttime=tr.stats.starttime,
                endtime=tr.stats.endtime,
            )
        else:
            my_log.warning(
                f"No station response info in stationXML file. Getting station response for {tr.id} from IRIS"
            )
            client = FDSNClient("IRIS")
            if check_FDSN(tr, client):
                sleep(0.25)
                inv = client.get_stations(
                    network=tr.stats.network,
                    station=tr.stats.station,
                    location=tr.stats.location,
                    channel=tr.stats.channel,
                    starttime=tr.stats.starttime,
                    endtime=tr.stats.endtime,
                    level="response",
                )
            else:
                my_log.warning(f"No data available for request. Removing {tr.id}")
                st.remove(tr)
                continue
        tr.stats.coordinates = inv.get_coordinates(tr.id, tr.stats.starttime)
        tr.inventory = inv

    return st


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


def grab_data(client, NSLC, T1, T2, fill_value=0):
    """
    Retrieve waveform data for specified channels and time range.

    Args:
        client (ObsPy Client): Client object to fetch data.
        NSLC (list or dict): List of channel identifiers (e.g., 'NET.STA.LOC.CHA').
        T1 (UTCDateTime): Start time for data retrieval.
        T2 (UTCDateTime): End time for data retrieval.
        fill_value (int or str): Value to fill gaps in data (default is 0).

    Returns:
        Stream: ObsPy Stream object containing the retrieved data.
    """
    my_log.info(f"Grabbing data from {client.name}...")

    st = Stream()

    if isinstance(NSLC, dict):
        NSLC = list(NSLC.keys())

    for nslc in NSLC:
        nslc = nslc.replace("--", "")  # Remove placeholder for empty location codes
        try:
            # Fetch waveform data for the specified channel and time range
            tr = client.get_waveforms(*nslc.split('.'), T1, T2)
            if len(tr) > 1:
                # Handle cases with multiple traces (e.g., due to gaps)
                if fill_value == 0 or fill_value is None:
                    tr.detrend("demean")
                    tr.taper(max_percentage=0.01)
                for sub_trace in tr:
                    # Ensure consistent data types and sampling rates
                    if sub_trace.data.dtype.name != "int32":
                        sub_trace.data = sub_trace.data.astype("int32")
                    if sub_trace.stats.sampling_rate != np.round(sub_trace.stats.sampling_rate):
                        sub_trace.stats.sampling_rate = np.round(sub_trace.stats.sampling_rate)
                my_log.info("Merging gappy data...")
                tr.merge(fill_value=fill_value)

            # Handle cases where the trace length is shorter than expected
            if tr[0].stats.endtime - tr[0].stats.starttime < T2 - T1:
                tr.detrend('demean')
                tr.taper(max_percentage=0.01)
        except Exception as e:
            my_log.error(f"Error occurred while grabbing data: {e}")
            my_log.warning(f"No data available for {nslc} from {client.name}. Creating empty Stream object.")
            tr = Stream()  # Create an empty stream if data retrieval fails

        # If no data is available, create a blank trace
        if not tr:
            from obspy import Trace
            from numpy import zeros
            tr = Trace()
            tr.id = nslc
            tr.stats['sampling_rate'] = 100
            tr.stats['starttime'] = T1
            tr.data = zeros(int((T2 - T1) * tr.stats["sampling_rate"]), dtype="int32")
        st += tr

    # Trim the stream to the specified time range and fill gaps
    st.trim(T1, T2, pad=True, fill_value=0)
    my_log.info("Detrending data...")
    st.detrend("demean")
    return st


def web_folders(t2, config, params):
    """
    Create directory structure for web output.

    Args:
        t2 (UTCDateTime): Current time.
        config (dict): Configuration dictionary.
        params (dict): Array parameters.

    Returns:
        None
    """
    # Create the base output directory if it doesn't exist
    out_web_dir = Path(config["OUT_WEB_DIR"])
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
        t2 (UTCDateTime): Current time.
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
    t1 = t2 - config[name]["DURATION"]
    tmp_name = name.replace(' ', '_')

    # Create output directories
    out_ascii_dir = Path(config["OUT_ASCII_DIR"])
    array_dir = out_ascii_dir / tmp_name
    month_dir = array_dir / t1.strftime('%Y-%m')
    month_dir.mkdir(parents=True, exist_ok=True)

    filename = month_dir / f"{tmp_name}_{t1.strftime('%Y-%m-%d')}.txt"

    # Adjust azimuth values to be within 0-360 degrees
    azimuth[azimuth < 0] += 360

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
        t2 (UTCDateTime): Current time.
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
    out_valve_dir = Path(config["OUT_VALVE_DIR"])
    out_valve_dir.mkdir(parents=True, exist_ok=True)
    filename = out_valve_dir / f"{name}_{t2.strftime('%Y%m%d-%H%M')}.txt"
    A.to_csv(filename, index=False, header=True, sep=',', float_format='%.3f')
    return

