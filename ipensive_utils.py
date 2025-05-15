import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from obspy import Stream, UTCDateTime, read_inventory
from obspy.clients.earthworm import Client as EWClient
from obspy.clients.fdsn import Client as FDSNClient
from obspy.clients.filesystem.sds import Client as SDSClient
from obspy.clients.seedlink import Client as SLClient
from obspy.geodetics.base import gps2dist_azimuth
from obspy.core.util import AttribDict
import yaml

####### plotting imports #######
import matplotlib as m
m.use('Agg')  # Use a non-interactive backend for matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib import dates
from copy import deepcopy
from collections import Counter
rcParams.update({'font.size': 10})  # Set default font size for plots
################################


def load_config(config_file):
    """
    Load configuration from a YAML file and initialize network and array settings.

    Args:
        config_file (str): Path to the configuration file.

    Returns:
        dict: Configuration dictionary with additional metadata.
    """
    with open(config_file, "r") as file:
        config = yaml.safe_load(file)

    # Extract network and array information
    all_nets = list(config["NETWORKS"].keys())
    array_list = []
    for net in all_nets:
        for array in config["NETWORKS"][net]:
            config[array]["NETWORK_NAME"] = net
            config[array]["ARRAY_NAME"] = array
            array_list.append(array)

    config["network_list"] = all_nets
    config["array_list"] = array_list

    # Load data source configuration
    with open(config["DATA_SOURCE"], "r") as file:
        data_source = yaml.safe_load(file)

    client = get_obspy_client(data_source)
    
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
    if config["CLIENT_TYPE"].lower() == "fdsn":
        client = FDSNClient(config["HOSTNAME"])
        client.name = config["HOSTNAME"]

    elif config["CLIENT_TYPE"].lower() == "local_fdsn":
        client = FDSNClient("IRIS", service_mappings={"dataselect": config["LOCAL_FDSN"]})
        client.name = config["LOCAL_FDSN"]

    elif config["CLIENT_TYPE"].lower() == "sds":
        client = SDSClient(config["DIRECTORY"])
        if "FMT" in list(config.keys()):
            client.FMTSTR = config["FMT"]
        client.name = config["DIRECTORY"]

    elif config["CLIENT_TYPE"].lower() == "earthworm":
        client = EWClient(config["HOSTNAME"], config["PORT"])
        client.name = config["HOSTNAME"]

    elif config["CLIENT_TYPE"].lower() == "seedlink":
        if "PORT" in list(config.keys()):
            client = SLClient(config["HOSTNAME"], config["PORT"])
        else:
            client = SLClient(config["HOSTNAME"])
        
    return client


def write_to_log(day, config):
    """
    Write logs to a file for a specific day.

    Args:
        day (str): Date string in UTC.
        config (dict): Configuration dictionary.
    """
    # Determine the logs directory
    if 'LOGS_DIR' in dir(config) and len(config.LOGS_DIR) > 0:
        logs_dir = Path(config.LOGS_DIR)
    else:
        logs_dir = Path(__file__).parent / 'logs'

    # Create year and month directories if they don't exist
    year = UTCDateTime(day).strftime('%Y')
    month = UTCDateTime(day).strftime('%Y-%m')
    year_dir = logs_dir / year
    month_dir = year_dir / month

    year_dir.mkdir(parents=True, exist_ok=True)
    month_dir.mkdir(parents=True, exist_ok=True)

    # Create the log file
    log_file = month_dir / f"{UTCDateTime(day).strftime('%Y-%m-%d')}.log"
    log_file.touch()
    with log_file.open("a") as f:
        sys.stdout = sys.stderr = f
    return



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
        inventory = client.get_stations(
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
            print(f"Skipping metadata for {tr.id} due to missing data")
            tr.stats.coordinates = AttribDict({
                'latitude': np.nan,
                'longitude': np.nan,
                'elevation': 0.0
            })
            continue
        print(f"Getting metadata for {tr.id}")
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
            print(
                f"No station response info in stationXML file. Getting station response for {tr.id} from IRIS"
            )
            client = FDSNClient("IRIS")
            if check_FDSN(tr, client):
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
                print(f"No data available for request. Removing {tr.id}")
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
    lon0 = np.mean([tr.stats.coordinates.longitude for tr in st])
    lat0 = np.mean([tr.stats.coordinates.latitude for tr in st])
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
    print(f"Grabbing data from {client.name}...")

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
                print("Merging gappy data...")
                tr.merge(fill_value=fill_value)

            # Handle cases where the trace length is shorter than expected
            if tr[0].stats.endtime - tr[0].stats.starttime < T2 - T1:
                tr.detrend('demean')
                tr.taper(max_percentage=0.01)
        except:
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
    print("Detrending data...")
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


def missing_elements_adjust(st, skip_chans, array):
    """
    Adjust the indices of an array to account for missing channels.

    Args:
        st (Stream): ObsPy Stream object containing traces.
        skip_chans (list): List of channels to skip.
        array (numpy.ndarray): Array of indices to adjust.

    Returns:
        tuple: Adjusted array and a list of missing indices.
    """
    skip_chans.sort()  # Ensure the skipped channels are sorted
    missing_inds = []  # List to store indices of missing channels
    for chan in skip_chans:
        # Find the index of the missing channel in the stream
        ind = np.where(np.array([tr.id for tr in st]) == chan)[0]
        # Increment indices in the array that are greater than the missing index
        array[array > ind] += 1
        # Add the missing index (adjusted by +1) to the list
        missing_inds.append(ind + 1)
    return array, missing_inds


def plot_results(t1, t2, t, st, mccm, velocity, azimuth, lts_dict, config, array_params, skip_chans):
    """
    Generate and save plots for the results of array processing.

    Args:
        t1 (UTCDateTime): Start time of the data window.
        t2 (UTCDateTime): End time of the data window.
        t (numpy.ndarray): Array of timestamps for the results.
        st (Stream): ObsPy Stream object containing traces.
        mccm (numpy.ndarray): MCCM values for the results.
        velocity (numpy.ndarray): Trace velocities for the results.
        azimuth (numpy.ndarray): Back-azimuth values for the results.
        lts_dict (dict): Dictionary of LTS (Least Trimmed Squares) results.
        config (dict): Configuration dictionary.
        array_params (dict): Parameters for the array being processed.
        skip_chans (list): List of channels to skip.

    Returns:
        None
    """
    # Define output directories for plots
    d0 = config["OUT_WEB_DIR"] + '/' + array_params["NETWORK_NAME"] + '/' + array_params["ARRAY_NAME"] + '/' + str(t2.year)
    d2 = d0 + '/' + '{:03d}'.format(t2.julday)

    # Generate a time vector for plotting waveforms
    tvec = np.linspace(
        dates.date2num(st[0].stats.starttime.datetime),
        dates.date2num(st[0].stats.endtime.datetime),
        len(st[0].data),
    )
    T1 = dates.date2num(t1.datetime)
    T2 = dates.date2num(t2.datetime)

    # Define colormap and color axis limits for MCCM values
    cm = "RdYlBu_r"
    cax = 0.2, 1  # Colorbar range for MCCM values

    # Determine the layout of the plot based on whether MCCM is being plotted
    if array_params["PLOT_MCCM"]:
        ax_list = [["wave"], ["cc"], ["vel"], ["baz"], ["stas"]]
        hr_list = [1, 1, 1, 1, 0.66]
    else:
        ax_list = [["wave"], ["vel"], ["baz"], ["stas"]]
        hr_list = [1, 1, 1, 0.66]

    # Set default plot size and styling parameters
    size = (8, 10.5)
    trace_lw = 0.6
    scatter_lw = 0.3
    s_dot = 36
    hline_lw = 1
    wm_font = 18
    
    # Loop through plot sizes (big and small) to generate both full-size and thumbnail plots
    for plot_size in ["big", "small"]:
        if plot_size == "small":
            # Adjust parameters for thumbnail plots
            size = (2.1, 2.75)
            trace_lw = 0.1
            scatter_lw = 0.1
            hline_lw = 0.4
            s_dot = 8
            wm_font = 8

        # Create the figure and axes using a mosaic layout
        fig, ax = plt.subplot_mosaic(
            ax_list,
            figsize=size,
            height_ratios=hr_list
        )

        ############ Plot Waveforms ##############
        ##########################################
        ax["wave"].set_title(array_params["ARRAY_NAME"] + " " + array_params["ARRAY_LABEL"] + " Array")
        ax["wave"].plot(tvec, st[0].data, "k", linewidth=trace_lw)
        ax["wave"].axis("tight")
        ax["wave"].set_xlim(T1, T2)
        ymax = np.abs(list(ax["wave"].get_ylim())).max()
        ax["wave"].set_ylim(-ymax, ymax)
        if plot_size == "big":
            ax["wave"].xaxis_date()
            ax["wave"].fmt_xdata = dates.DateFormatter("%HH:%MM")
            ax["wave"].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
            ax["wave"].set_xticklabels([])
            ax["wave"].tick_params(direction="in", axis="x", top="on")
            ax["wave"].set_ylabel("Pressure [Pa]")
            ax["wave2"] = ax["wave"].twinx()
            ax["wave2"].set_yticks([])
            ax["wave2"].set_ylabel(
                f"{array_params['FREQMIN']:.1f} - {array_params['FREQMAX']:.1f} Hz",
                labelpad=6,
            )
        else:
            ax["wave"].set_xticks([])
            ax["wave"].set_yticks([])
        ##########################################


        ############ Plot cc values ##############
        ##########################################
        # Check if MCCM (Mean Cross-Correlation Metric) plotting is enabled
        if array_params["PLOT_MCCM"]:
            # Scatter plot of MCCM values with color mapping
            sc = ax["cc"].scatter(t, mccm, c=mccm, s=s_dot, edgecolors="k", lw=scatter_lw, cmap=cm)
            # Add a horizontal line for the MCCM threshold
            ax["cc"].axhline(array_params["MCTHRESH"], ls="--", lw=hline_lw, color="gray")
            # Adjust axis limits and scaling
            ax["cc"].axis("tight")
            ax["cc"].set_xlim(T1, T2)
            ax["cc"].set_ylim(cax)
            sc.set_clim(cax)
            # Configure axis labels and formatting for larger plots
            if plot_size == "big":
                ax["cc"].xaxis_date()
                ax["cc"].fmt_xdata = dates.DateFormatter("%HH:%MM")
                ax["cc"].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
                ax["cc"].set_xticklabels([])
                ax["cc"].tick_params(direction="in", axis="x", top="on")
                ax["cc"].set_ylabel(r"$M_{d}CCM$")
            else:
                # Remove ticks and labels for smaller plots
                ax["cc"].set_xticks([])
                ax["cc"].set_yticks([])
        ##########################################

        ########## Plot Trace Velocities #########
        ##########################################
        # Highlight the velocity range of interest with a shaded region
        ax["vel"].axhspan(
            array_params["VEL_MIN"],
            array_params["VEL_MAX"],
            facecolor="gray",
            alpha=0.25,
            edgecolor=None,
        )

        # Scatter plot of trace velocities with color mapping based on MCCM values
        sc = ax["vel"].scatter(t, velocity, c=mccm, s=s_dot, edgecolors="k", lw=scatter_lw, cmap=cm)
        # Set y-axis limits based on the array type
        if array_params["ARRAY_LABEL"] == "Hydroacoustic":
            ax["vel"].set_ylim(1.2, 1.8)  # Typical range for hydroacoustic arrays
        else:
            ax["vel"].set_ylim(0.15, 0.6)  # Typical range for other arrays
        # Set x-axis limits to match the time range
        ax["vel"].set_xlim(T1, T2)
        # Set the color limits for the scatter plot
        sc.set_clim(cax)
        # Configure axis labels and formatting for larger plots
        if plot_size == "big":
            ax["vel"].xaxis_date()
            ax["vel"].fmt_xdata = dates.DateFormatter("%HH:%MM")
            ax["vel"].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
            ax["vel"].set_xticklabels([])
            ax["vel"].tick_params(direction="in", axis="x", top="on")
            ax["vel"].set_ylabel("Trace Velocity\n [km/s]")
        else:
            # Remove ticks and labels for smaller plots
            ax["vel"].set_xticks([])
            ax["vel"].set_yticks([])
        ##########################################

        ########### Plot Back-azimuths ###########
        ##########################################
        # Configure the style for the text boxes
        box_style = {'facecolor': 'white', 'edgecolor': 'white', 'pad': 0}
        # Adjust azimuth range if AZ_MAX is less than AZ_MIN
        az_min = deepcopy(array_params['AZ_MIN'])
        az_max = deepcopy(array_params['AZ_MAX'])
        tmp_azimuth = deepcopy(azimuth)
        if array_params['AZ_MAX'] < array_params['AZ_MIN']:
            az_min = az_min - 360
            # Loop through targets and adjust back-azimuths
            for target in array_params['TARGETS']:
                if isinstance(target, dict):
                    target = list(target.keys())[0]
                baz = array_params[target]
                if baz > 180:
                    # Plot horizontal lines and labels for targets
                    ax["baz"].axhline(baz - 360, ls='--', lw=hline_lw, color='gray', zorder=-1)
                    if plot_size == "big":
                        ax["baz"].text(t[1], baz - 360, target, bbox=box_style, fontsize=8, va='center', style='italic', zorder=10)
                else:
                    ax["baz"].axhline(baz, ls='--', lw=hline_lw, color='gray', zorder=-1)
                    if plot_size == "big":
                        ax["baz"].text(t[1], baz, target, bbox=box_style, fontsize=8, va='center', style='italic', zorder=10)
            tmp_azimuth[tmp_azimuth > 180] += -360
        else:
            # Plot horizontal lines and labels for targets
            for target in array_params['TARGETS']:
                baz = array_params[target]
                ax["baz"].axhline(baz, ls='--', lw=hline_lw, color='gray', zorder=-1)
                if plot_size == "big":
                    ax["baz"].text(t[1], baz, target, bbox=box_style, fontsize=8, va='center', style='italic', zorder=10)

        # Scatter plot for back-azimuth values
        sc = ax["baz"].scatter(t, tmp_azimuth, c=mccm, s=s_dot, edgecolors='k', lw=scatter_lw, cmap=cm, zorder=1000)
        ax["baz"].set_ylim(az_min, az_max)
        ax["baz"].set_xlim(T1, T2)
        sc.set_clim(cax)
        if plot_size == "big":
            # Configure x-axis for larger plots
            ax["baz"].xaxis_date()
            ax["baz"].fmt_xdata = dates.DateFormatter("%HH:%MM")
            ax["baz"].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
            ax["baz"].set_xticklabels([])
            ax["baz"].tick_params(direction="in", axis="x", top="on")
            ax["baz"].set_ylabel("Trace Velocity\n [km/s]")
        else:
            # Remove ticks for smaller plots
            ax["baz"].set_xticks([])
            ax["baz"].set_yticks([])
        ##########################################

        ########## Plot Dropped Channels #########
        ##########################################
        if len(lts_dict) == 0:
            # Handle case where there are not enough channels for LTS
            txt_str = "Not enough channels for LTS"
            ax["stas"].text(0.5, 0.5, txt_str, transform=ax["stas"].transAxes, color='grey', alpha=0.7, fontsize=wm_font, va="center", ha="center")
            n = len(st)
            missing_inds = []
            # Identify missing channels
            for chan in skip_chans:
                ind = np.where(np.array([tr.id for tr in st]) == chan)[0] + 1
                missing_inds.append(ind)
        else:
            # Process LTS dictionary to plot dropped channels
            ndict = deepcopy(lts_dict)
            n = ndict['size']
            ndict.pop('size', None)
            tstamps = list(ndict.keys())
            tstampsfloat = [float(ii) for ii in tstamps]
            cm2 = plt.get_cmap('binary', (n - 1))
            for jj in range(len(tstamps)):
                z = Counter(list(ndict[tstamps[jj]]))
                keys, vals = z.keys(), z.values()
                keys, vals = np.array(list(keys)), np.array(list(vals))
                pts = np.tile(tstampsfloat[jj], len(keys))
                keys, missing_inds = missing_elements_adjust(st, skip_chans, keys)
                # Scatter plot for dropped channels
                sc_stas = ax["stas"].scatter(
                    pts,
                    keys,
                    c=vals,
                    s=s_dot,
                    edgecolors="k",
                    lw=scatter_lw,
                    cmap=cm2,
                    vmin=0.5,
                    vmax=n - 0.5,
                )
        # Configure y-axis for station labels
        ax["stas"].set_ylim(0, len(st) + 1)
        ax["stas"].set_yticks(np.arange(1, len(st) + 1))
        ax["stas"].set_yticklabels([f"{tr.stats.station}.{tr.stats.location}" for tr in st], fontsize=8)
        ax["stas"].set_xlim(T1, T2)

        if plot_size == "big":
            # Add scatter plot for missing channels
            if n > 3 and ndict:
                for m_ind in missing_inds:
                    ax["stas"].scatter(t, m_ind * np.ones(t.shape[0]), marker=".", color="indianred", s=10)
            # Configure x-axis for larger plots
            ax["stas"].xaxis_date()
            ax["stas"].tick_params(axis='x', labelbottom='on')
            ax["stas"].fmt_xdata = dates.DateFormatter('%HH:%MM')
            ax["stas"].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
            ax["stas"].set_xlabel('UTC Time [' + t1.strftime('%Y-%b-%d') + ']')
        else:
            # Remove ticks for smaller plots
            ax["stas"].set_xticks([])
            ax["stas"].set_yticks([])
        ax["stas"].invert_yaxis()
        ##########################################

        ########## Adjust & Save Figure ##########
        ##########################################
        if plot_size == "big":
            # Adjust layout and add colorbars for larger plots
            plt.subplots_adjust(left=0.1, right=0.9, top=0.97, bottom=0.05, hspace=0.1)

            ax_str = "cc" if array_params["PLOT_MCCM"] else "vel"
            ctop = ax[ax_str].get_position().y1
            cbot = ax["baz"].get_position().y0
            cbaxes_mccm = fig.add_axes([0.91, cbot, 0.02, ctop - cbot])
            hc = plt.colorbar(sc, cax=cbaxes_mccm)
            hc.set_label(r'$M_{d}CCM$')
            if n > 3 and array_params["LTS_ALPHA"] < 1:
                ctop = ax["stas"].get_position().y1
                cbot = ax["stas"].get_position().y0
                if ndict:
                    cbaxes_stas = fig.add_axes([0.91, cbot, 0.02, ctop - cbot])
                    hc_stas = plt.colorbar(sc_stas, cax=cbaxes_stas)
                    hc_stas.set_label(r"# Dropped Pairs")
                    hc_stas.set_ticks(np.arange(1, n))
                else:
                    txt_str = "No dropped channels"
                    ax["stas"].text(0.5, 0.5, txt_str, transform=ax["stas"].transAxes, color='grey', alpha=0.7, fontsize=wm_font, va="center", ha="center")
            
            # Save the full-size plot
            filename = d2 + '/' + array_params["ARRAY_NAME"] + '_' + t2.strftime('%Y%m%d-%H%M') + '.png'
            fig.savefig(filename, dpi=72, format='png')
            plt.close("all")
        else:
            # Adjust layout and save thumbnail plot
            plt.subplots_adjust(left=0, right=1, top=0.99, bottom=0.01, hspace=0.03)
            filename = d2 + '/' + array_params["ARRAY_NAME"] + '_' + t2.strftime('%Y%m%d-%H%M') + '_thumb.png'
            fig.savefig(filename, format='png', pad_inches=0, dpi=72)
            plt.close('all')
