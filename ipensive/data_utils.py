import logging
import numpy as np
from obspy import Stream
from obspy.clients.earthworm import Client as EWClient
from obspy.clients.fdsn import Client as FDSNClient
from obspy.clients.filesystem.sds import Client as SDSClient
from obspy.clients.seedlink import Client as SLClient


my_log = logging.getLogger(__name__)


def get_obspy_client(config):
    """
    Initialize an ObsPy client based on the configuration.

    Args:
        config (dict): Configuration dictionary for the client.

    Returns:
        ObsPy client object.
    """

    if "TIMEOUT" not in config:
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

    else:
        client = None
        my_log.error(f"CLIENT_TYPE {config['CLIENT_TYPE']} not recognized. Exiting.")

    return client


def grab_data(client, NSLC, T1, T2, fill_value=0):
    """
    Retrieve waveform data for specified channels and time range.

    Args:
        client (ObsPy Client): Client object to fetch data.
        NSLC (list or dict): List of channel identifiers (e.g., 'NET.STA.LOC.CHA').
        T1 (obspy.UTCDateTime): Start time for data retrieval.
        T2 (obspy.UTCDateTime): End time for data retrieval.
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


def preprocess_data(st, t1, t2, skip_chans, array_params):
    for tr in st:
        if tr.id in skip_chans:
            continue
        if isinstance(array_params["NSLC"], dict):
            tr.data = tr.data / array_params["NSLC"][tr.id]["gain"]
        else:
            tr.remove_sensitivity(tr.inventory)
    st.detrend("demean")
    st.taper(max_percentage=None, max_length=array_params["TAPER"])
    st.filter(
        "bandpass",
        freqmin=array_params["FREQMIN"],
        freqmax=array_params["FREQMAX"],
        corners=2,
        zerophase=True,
    )
    st.trim(t1, t2 + array_params["WINDOW_LENGTH"])


def QC_data(st, array_params):

    #### Check for enough data ####
    check_st = st.copy()
    skip_chans = []
    good_data = True
    for tr in check_st:
        if np.sum(np.abs(tr.data)) == 0:  # Check for blank traces
            skip_chans.append(tr.id)
            check_st.remove(tr)
    if len(check_st) < array_params["MIN_CHAN"]:
        my_log.warning("Too many blank traces. Skipping.")
        good_data = False
        return good_data, skip_chans
    ########################

    #### Check for gappy data ####
    for tr in check_st:
        if np.any([np.any(tr.data == 0)]):  # Check for gaps in data
            check_st.remove(tr)
    if len(check_st) < array_params["MIN_CHAN"]:
        my_log.warning("Too gappy. Skipping.")
        good_data = False

    return good_data, skip_chans
