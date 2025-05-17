from obspy.core.inventory.inventory import Inventory
from obspy.clients.fdsn import Client
from obspy import UTCDateTime as utc
from .ipensive_utils import write_to_log
import os
from time import sleep

import logging
my_log = logging.getLogger(__name__)

def get_stations(config):

    NSLC = []
    for net in list(config["NETWORKS"].keys()):
        for array in config["NETWORKS"][net]:
            NSLC += config[array]["NSLC"]

    return NSLC


def update_stationXML(config):

    if os.getenv("FROMCRON") == "yep":
        write_to_log(utc.utcnow(), config)
    
    client_iris = Client("IRIS")
    NSLC = get_stations(config)
    print("______ Begin Updating Metadata ______")
    print("______ " + utc.utcnow().strftime("%Y-%m-%d %H:%M:%S") + " ______")

    inventory = Inventory()
    for nslc in NSLC:

        print(nslc)
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
                print(f"Error on attempt number {attempts:g}:")
                print(f"\t{ex}")

    inventory.write(config["STATION_XML"], format="STATIONXML")

    print("^^^^^^ Finished Updating Metadata ^^^^^^")
    return