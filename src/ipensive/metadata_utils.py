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

    client_iris = FDSNClient("IRIS")
    NSLC = get_stations(config)
    my_log.info("______ Begin Updating Metadata ______")
    my_log.info("______ " + utc.utcnow().strftime("%Y-%m-%d %H:%M:%S") + " ______")

    inventory = Inventory()
    for nslc in NSLC:
        sleep(0.25)
        my_log.info(f"Updating metadata for {nslc}")
        net, sta, loc, chan = nslc.split(".")
        client = client_iris
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
    except Exception as err:
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


def add_metadata(st, config, array_name, skip_chans=[]):
    """
    Add metadata to traces in a stream.

    Args:
        st (Stream): ObsPy Stream object.
        config (dict): Configuration dictionary.
        skip_chans (list): List of channels (NSLC) to skip.

    Returns:
        Stream: Stream with updated metadata.
    """
    
    import warnings

    warnings.simplefilter("ignore", UserWarning, append=True)

    lat_list = []
    lon_list = []

    if isinstance(config[array_name]["NSLC"], dict):
        my_log.info(f"Adding coordinate info for {array_name} directly from config file")
        st = add_coordinate_info(st, config, array_name)
        for tr in st:
            lat_list.append(tr.stats.coordinates.latitude)
            lon_list.append(tr.stats.coordinates.longitude)
        return st, lat_list, lon_list

    if "STATION_XML" in config.keys():
        my_log.info(f"Adding metadata from {config['STATION_XML']}")
        inventory = read_inventory(config["STATION_XML"])

    empty_coords = AttribDict({
                    'latitude': np.nan,
                    'longitude': np.nan,
                    'elevation': np.nan
                })

    for tr in st:
        my_log.info(f"Getting metadata for {tr.id}")

        if tr.id in skip_chans:
            my_log.info(f"{tr.id} is in the skip list. Adding empty coordinates.")
            tr.stats.coordinates = empty_coords

        elif check_inventory(tr, inventory):
            inv = inventory.select(
                network=tr.stats.network,
                station=tr.stats.station,
                location=tr.stats.location,
                channel=tr.stats.channel,
                starttime=tr.stats.starttime,
                endtime=tr.stats.endtime,
            )
            tr.stats.coordinates = inv.get_coordinates(tr.id, tr.stats.starttime)
            tr.inventory = inv

        else:
            my_log.warning(
                f"No station response info in stationXML file. Getting station response for {tr.id} from IRIS"
            )

            client = FDSN_connect("IRIS", max_tries=3)
            if not client:
                my_log.error(f"IRIS FDSN client unavailable for channel {tr.id}")
                my_log.warning("...Adding empty coordinates. This might break things")
                tr.stats.coordinates = empty_coords
            elif check_FDSN(tr, client):
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
                tr.stats.coordinates = inv.get_coordinates(tr.id, tr.stats.starttime)
                tr.inventory = inv
            else:
                my_log.error(f"No data available for request for channel {tr.id}. Adding NaNs")
                my_log.warning("...This might break things")
                tr.stats.coordinates = empty_coords

        lat_list.append(tr.stats.coordinates.latitude)
        lon_list.append(tr.stats.coordinates.longitude)

    return st, lat_list, lon_list