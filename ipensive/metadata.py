import logging
from time import sleep
from obspy.core.inventory.inventory import Inventory
from obspy.clients.fdsn import Client
from obspy import UTCDateTime as utc


my_log = logging.getLogger(__name__)

def get_stations(config):

    NSLC = []
    for net in list(config["NETWORKS"].keys()):
        for array in config["NETWORKS"][net]:
            NSLC += config[array]["NSLC"]

    return NSLC


def update_stationXML(config):

    client_iris = Client("IRIS")
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