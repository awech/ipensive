import logging
from time import sleep
import numpy as np
from obspy.core.inventory.inventory import Inventory
from obspy import UTCDateTime as utc
from obspy import read_inventory
from obspy.core.util import AttribDict
from obspy.clients.fdsn import Client as FDSNClient


my_log = logging.getLogger(__name__)


def get_stations(config):
    """Get a list of station NSLC codes from the configuration.

    Args:
        config (dict): Configuration dictionary.

    Returns:
        list: List of station NSLC codes.
    """

    NSLC = []
    for array in config["array_list"]:
        NSLC += config[array]["NSLC"]

    return NSLC


def update_stationXML(config):
    """Update the station metadata XML file.

    Args:
        config (dict): Configuration dictionary.
    """

    client_earthscope = FDSNClient("earthscope")
    NSLC = get_stations(config)
    my_log.info("______ Begin Updating Metadata ______")
    my_log.info("______ " + utc.utcnow().strftime("%Y-%m-%d %H:%M:%S") + " ______")

    inventory = Inventory()
    for nslc in NSLC:
        sleep(0.25)
        my_log.info(f"Updating metadata for {nslc}")
        net, sta, loc, chan = nslc.split(".")
        client = client_earthscope
        attempts = 0
        while attempts < 4:
            try:
                inventory += client.get_stations(
                    station=sta,
                    network=net,
                    channel=chan,
                    location=loc,
                    level="response",
                    format="xml",
                )
                break
            except Exception as ex:
                sleep(1)
                attempts += 1
                my_log.warning(f"Error on attempt number {attempts:g}:")
                my_log.error(f"\t{ex}")

    inventory.write(config["STATION_XML"], format="STATIONXML")

    my_log.info("^^^^^^ Finished Updating Metadata ^^^^^^\n")
    return


def FDSN_connect(client_name, max_tries=3):

    client = []
    attempts = 0
    while attempts < max_tries:
        try:
            client = FDSNClient(client_name, timeout=10)
            break
        except Exception as ex:
            sleep(1)
            attempts += 1
            my_log.warning(f"Error on attempt number {attempts:g}")
            my_log.error(ex)

    return client


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
    except Exception as err:  # pragma: no cover
        if "No data available for request." in err.args[0]:
            value = False
    return value


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


def get_inventory(tr, inventory, skip_chans=[]):
    """
    Get the subset inventory (with response) for a single trace.

    Looks up the trace's channel/epoch in the provided (already loaded)
    inventory. Falls back to an Earthscope FDSN request if the channel
    isn't found locally. Channels in ``skip_chans`` are skipped entirely
    (no local lookup, no network fallback) since they're already known
    to be unusable (e.g. dead/blank channels).

    Args:
        tr (Trace): ObsPy Trace object.
        inventory (Inventory): Pre-loaded station inventory (e.g. from STATION_XML).
        skip_chans (list): List of channels (NSLC) to skip.

    Returns:
        Inventory or None: The subset inventory for this trace's channel/epoch,
            or None if the channel should be skipped or no response info could
            be found.
    """

    if tr.id in skip_chans:
        my_log.info(f"{tr.id} is in the skip list. Skipping inventory lookup.")
        return None

    if check_inventory(tr, inventory):
        return inventory.select(
            network=tr.stats.network,
            station=tr.stats.station,
            location=tr.stats.location,
            channel=tr.stats.channel,
            starttime=tr.stats.starttime,
            endtime=tr.stats.endtime,
        )

    my_log.warning( # pragma: no cover
        f"No station response info in stationXML file. Getting station response for {tr.id} from Earthscope"
    )

    client = FDSN_connect("earthscope", max_tries=3) # pragma: no cover
    if not client: # pragma: no cover
        my_log.error("Earthscope FDSN client unavailable")
        return None
    elif check_FDSN(tr, client): # pragma: no cover
        sleep(0.25)
        return client.get_stations(
            network=tr.stats.network,
            station=tr.stats.station,
            location=tr.stats.location,
            channel=tr.stats.channel,
            starttime=tr.stats.starttime,
            endtime=tr.stats.endtime,
            level="response",
        )
    else: # pragma: no cover
        my_log.error(f"No data available for request for channel {tr.id}.")
        return None


def add_coordinates_from_config(st, config, array_name):
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


def add_coordinates(st, config, array_name, skip_chans=[]):
    """
    Add coordinate (and inventory) metadata to traces in a stream.

    For each trace, looks up its station coordinates and (when using a
    STATION_XML-based array) stashes the resolved inventory subset on
    ``tr.inventory`` for later use (e.g. by ``remove_gain``), so the
    inventory lookup/Earthscope fallback only happens once per trace.

    Args:
        st (Stream): ObsPy Stream object.
        config (dict): Configuration dictionary.
        array_name (str): Name of the array.
        skip_chans (list): List of channels (NSLC) to skip.

    Returns:
        tuple: (st, lat_list, lon_list)
    """

    import warnings

    warnings.simplefilter("ignore", UserWarning, append=True)

    lat_list = []
    lon_list = []

    if isinstance(config[array_name]["NSLC"], dict):
        my_log.info(f"Adding coordinate info for {array_name} directly from config file")
        st = add_coordinates_from_config(st, config, array_name)
        for tr in st:
            lat_list.append(tr.stats.coordinates.latitude)
            lon_list.append(tr.stats.coordinates.longitude)
        return st, lat_list, lon_list

    if "STATION_XML" in config[array_name].keys():
        my_log.info(f"Adding metadata from {config[array_name]['STATION_XML']}")
        inventory = read_inventory(config[array_name]["STATION_XML"])

    empty_coords = AttribDict({
                    'latitude': np.nan,
                    'longitude': np.nan,
                    'elevation': np.nan
                })

    for tr in st:
        my_log.info(f"Getting metadata for {tr.id}")

        inv = get_inventory(tr, inventory, skip_chans)
        if inv is None:
            if tr.id not in skip_chans: # pragma: no cover
                my_log.warning("...Adding empty coordinates. This might break things")
            tr.stats.coordinates = empty_coords
        else:
            tr.stats.coordinates = inv.get_coordinates(tr.id, tr.stats.starttime)
            tr.stats.inventory = inv

        lat_list.append(tr.stats.coordinates.latitude)
        lon_list.append(tr.stats.coordinates.longitude)

    return st, lat_list, lon_list


def remove_gain(st, array_params):
    """
    Remove instrument gain/sensitivity from traces in a stream.

    Uses the manual per-channel gain value for manually-configured (dict)
    NSLC arrays, or the inventory attached to each trace (via
    ``add_coordinates``) for STATION_XML-based arrays. Traces with no
    attached inventory (e.g. skipped/dead channels) are left as-is.

    Note: inventory is stashed on ``tr.stats.inventory`` (not a bare
    ``tr.inventory`` attribute) since ``Trace.stats`` is deep-copied by
    ``Stream.merge()`` and ``Trace.copy()``, while arbitrary attributes
    set directly on a ``Trace`` object are not preserved by either.

    Args:
        st (Stream): ObsPy Stream object.
        array_params (dict): Array parameters.

    Returns:
        Stream: Stream with gain removed.
    """
    for tr in st:
        if isinstance(array_params["NSLC"], dict):
            tr.data = tr.data / array_params["NSLC"][tr.id]["gain"]
        elif "inventory" in tr.stats and tr.stats.inventory is not None:
            tr.remove_sensitivity(tr.stats.inventory)
        else: # pragma: no cover
            my_log.warning(f"{tr.id}: no inventory attached. Skipping gain removal.")
    return st