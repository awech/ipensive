import logging
import numpy as np
from matplotlib import dates
from obspy import Stream
from obspy import UTCDateTime as utc
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

    if "TIMEOUT" not in config: # pragma: no cover
        config["TIMEOUT"] = 30

    if config["CLIENT_TYPE"].lower() == "fdsn":
        client = FDSNClient(config["HOSTNAME"], timeout=config["TIMEOUT"])
        client.name = config["HOSTNAME"]

    elif config["CLIENT_TYPE"].lower() == "local_fdsn": # pragma: no cover
        client = FDSNClient("earthscope", service_mappings={"dataselect": config["LOCAL_FDSN"]}, timeout=config["TIMEOUT"])
        client.name = config["LOCAL_FDSN"]

    elif config["CLIENT_TYPE"].lower() == "sds": # pragma: no cover
        client = SDSClient(config["DIRECTORY"])
        if "FMT" in list(config.keys()):
            client.FMTSTR = config["FMT"]
        client.name = config["DIRECTORY"]

    elif config["CLIENT_TYPE"].lower() == "earthworm": # pragma: no cover
        client = EWClient(config["HOSTNAME"], config["PORT"], timeout=config["TIMEOUT"])
        client.name = config["HOSTNAME"]

    elif config["CLIENT_TYPE"].lower() == "seedlink": # pragma: no cover
        if "PORT" in list(config.keys()):
            client = SLClient(config["HOSTNAME"], config["PORT"], timeout=config["TIMEOUT"])
        else:
            client = SLClient(config["HOSTNAME"], timeout=config["TIMEOUT"])
        client.name = config["HOSTNAME"]

    else: # pragma: no cover
        client = None
        my_log.error(f"CLIENT_TYPE {config['CLIENT_TYPE']} not recognized. Exiting.")

    return client


def grab_data(client, NSLC, T1, T2):
    """
    Retrieve waveform data for specified channels and time range.

    Args:
        client (ObsPy Client): Client object to fetch data.
        NSLC (list or dict): List of channel identifiers (e.g., 'NET.STA.LOC.CHA').
        T1 (obspy.UTCDateTime): Start time for data retrieval.
        T2 (obspy.UTCDateTime): End time for data retrieval.

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
            if len(tr) > 1: # pragma: no cover
                # Handle cases with multiple traces (e.g., due to gaps)
                for sub_trace in tr:
                    # Ensure consistent data types and sampling rates
                    sub_trace = _qc_sub_trace(sub_trace)
                if not tr.get_gaps():
                    # handle case where multiple traces returned with no gaps between them
                    my_log.info(f"{nslc}: Multiple traces returned with no gaps between. Simple merge")
                    tr.merge()

        except Exception as e: # pragma: no cover
            my_log.error(f"Error occurred while grabbing data: {e}")
            my_log.warning(f"No data available for {nslc} from {client.name}. Creating empty Stream object.")
            tr = Stream()  # Create an empty stream if data retrieval fails

        # If no data is available, create a blank trace
        if not tr: # pragma: no cover
            from obspy import Trace
            from numpy import zeros
            tr = Trace()
            tr.id = nslc
            tr.stats['sampling_rate'] = 100
            tr.stats['starttime'] = T1
            tr.data = zeros(int((T2 - T1) * tr.stats["sampling_rate"]), dtype="int32")
        st += tr

    return st


def preprocess_data(ST, t1, t2, array_params):
    """
    Preprocess seismic data by tapering, filtering, gap handling, and trimming.

    Args:
        ST (obspy.Stream): Stream containing seismic traces.
        t1 (obspy.UTCDateTime): Start time for trimming.
        t2 (obspy.UTCDateTime): End time for trimming.
        array_params (dict): Array parameters including filtering and tapering settings.

    Returns:
        obspy.Stream: Preprocessed stream.
    """

    st = ST.copy()

    st.detrend("demean")
    st.taper(max_percentage=None, max_length=array_params["TAPER"])
    st.filter(
        "bandpass",
        freqmin=array_params["FREQMIN"],
        freqmax=array_params["FREQMAX"],
        corners=2,
        zerophase=True,
    )
    gaps = st.get_gaps()
    if gaps:
        my_log.warning(f"Gappy data: {len(gaps)} gap(s)/overlap(s)")
        for net, sta, loc, chan, t_last, t_next, delta, samples in gaps:
            my_log.warning(
                f"{net}.{sta}.{loc}.{chan}: {t_last} -> {t_next} "
                f"(delta={delta:.3f}s, samples={samples})"
            )
        my_log.warning("Attempting to merge (fill_value=0)")
        st.merge(fill_value=0)
        
    st.trim(t1, t2 + array_params["WINDOW_LENGTH"], pad=True, fill_value=0)

    return st


def _qc_sub_trace(sub_trace):
    """
    Ensure consistent data type and sampling rate for a trace.

    Casts trace data to int32 if not already, and rounds the sampling
    rate to the nearest integer if it is not already a whole number.

    Args:
        sub_trace (obspy.Trace): A single seismic trace to check.

    Returns:
        obspy.Trace: The trace with corrected data type and sampling rate.
    """
    if sub_trace.data.dtype.name != "int32":
        my_log.info(f"{sub_trace.id}: changing dtype to int32")
        sub_trace.data = sub_trace.data.astype("int32")
    if sub_trace.stats.sampling_rate != np.round(sub_trace.stats.sampling_rate):
        my_log.info(f"{sub_trace.id}: sampling rate is non-integer. Rounding to nearest integer value")
        sub_trace.stats.sampling_rate = np.round(sub_trace.stats.sampling_rate)
    return sub_trace


def QC_data(st, array_params):
    """
    Quality control for seismic data.

    Blank channels or channels with a fraction of zero-filled (gap) 
    samples exceeding MAX_GAP_FRACTION are flagged and skipped.

    Args:
        st (obspy.Stream): Stream containing seismic traces.
        array_params (dict): Array parameters including quality control thresholds.

    Returns:
        tuple: (good_data, skip_chans) where good_data is a boolean indicating
                if the data passed QC and skip_chans is a list of channels to skip.
    """

    max_gap_fraction = array_params.get("MAX_GAP_FRACTION", 0.5)

    #### Check for enough data ####
    check_st = st.copy()
    skip_chans = []
    good_data = True
    for tr in check_st:
        if np.sum(np.abs(tr.data)) == 0: # pragma: no cover
            # Check for blank traces
            skip_chans.append(tr.id)
            check_st.remove(tr)
    if len(check_st) < array_params["MIN_CHAN"]: # pragma: no cover
        my_log.warning("Too many blank traces. Skipping.")
        good_data = False
        return good_data, skip_chans
    ########################

    #### Check for gappy data ####
    for tr in check_st:
        gap_fraction = np.count_nonzero(tr.data == 0) / tr.stats.npts
        if gap_fraction > max_gap_fraction: # pragma: no cover
            # Gap exceeds tolerance. Flag channel for exclusion.
            my_log.warning(f"{tr.id}: {gap_fraction:.1%} gap exceeds MAX_GAP_FRACTION ({max_gap_fraction:.1%}). Skipping channel.")
            skip_chans.append(tr.id)
            check_st.remove(tr)
    if len(check_st) < array_params["MIN_CHAN"]: # pragma: no cover
        my_log.warning("Too gappy. Skipping.")
        good_data = False

    return good_data, skip_chans


def get_pressures(st, t, array_params, skip_chans=[]):
    """Extract pressure data from the seismic stream.

    Uses the median peak amplitude across all channels not in skip_chans,
    so a single dead/excluded channel can't bias the reported pressure
    (either by being picked directly, or by skewing a mean). Note:
    PLOTCHAN (which selects a single channel for waveform plotting in
    plotting_utils.py) intentionally does NOT affect this calculation -
    pressure is a whole-array data output, not a plotting concern.

    Args:
        st (obspy.Stream): Stream containing seismic traces.
        t (np.ndarray): Array of time values (matplotlib dates) from ltsva.
        array_params (dict): Array parameters including window length.
        skip_chans (list): List of channels (NSLC) excluded by QC (e.g. dead
            or too-gappy channels) that should not contribute to the pressure
            estimate.

    Returns:
        np.ndarray: Array of pressure values.
    """

    good_st = Stream([tr for tr in st if tr.id not in skip_chans])
    keep_st = good_st if len(good_st) > 0 else st  # pragma: no cover

    pressure = []
    for ti in t:
        t1 = utc(dates.num2date(ti)) - array_params["WINDOW_LENGTH"] / 2
        t2 = t1 + array_params["WINDOW_LENGTH"]
        peak_amps = [np.max(np.abs(tr.slice(t1, t2).data)) for tr in keep_st]
        pressure.append(np.median(peak_amps))
    pressure = np.array(pressure)
    return pressure