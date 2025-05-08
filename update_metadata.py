from obspy.core.inventory.inventory import Inventory
from obspy.clients.fdsn import Client
from obspy import UTCDateTime as utc
from ipensive_utils import write_to_log
import os
import argparse
from time import sleep
import yaml


def get_stations(config):

    SCNL = []
    for net in list(config["NETWORKS"].keys()):
        for array in config["NETWORKS"][net]:
            SCNL += config[array]["SCNL"]

    return SCNL


def update_stationXML(config):

    if os.getenv("FROMCRON") == "yep":
        write_to_log(utc.utcnow(), config)

    client_iris = Client("IRIS")
    SCNL = get_stations(config)
    print("______ Begin Updating Metadata ______")
    print("______ " + utc.utcnow().strftime("%Y-%m-%d %H:%M:%S") + " ______")

    inventory = Inventory()
    for scnl in SCNL:

        print(scnl)
        sta, chan, net, loc = scnl.split(".")
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        epilog="e.g.: python update_metadata.py -c <filename.yml>"
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        help="Name of the config file (yml)",
        default="config.yml",
    )
    args = parser.parse_args()
    config_file = args.config
    with open(config_file, "r") as file:
        config = yaml.safe_load(file)
    update_stationXML(config)
