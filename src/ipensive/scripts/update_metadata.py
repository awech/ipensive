import argparse
import numpy as np
from obspy import UTCDateTime as utc
import ipensive.ipensive_utils as utils
from ipensive.metadata_utils import update_stationXML

def main():
    """
    Main entry point for the metadata update script.
    """

    parser = argparse.ArgumentParser(
        epilog="e.g.: python update_metadata.py -c <filename.yml>"
    )

    parser.add_argument(
        "-c",
        "--config",
        type=str,
        help="Name of the config file (yml)",
        default=utils.get_config_file(),
    )
    args = parser.parse_args()
    config_file = args.config
    config = utils.load_ipensive_config(config_file)

    utils.setup_logging(utc.utcnow(), config)

    xml_files = [config[array]["STATION_XML"] for array in config["array_list"]]
    for tmp_xml in np.unique(xml_files):
        tmp_config = config.copy()
        tmp_array_list = []
        for array in tmp_config["array_list"]:
            if tmp_config[array]["STATION_XML"] != tmp_xml:
                tmp_config.pop(array)
            else:
                tmp_array_list.append(array)
        tmp_config["array_list"] = tmp_array_list
        tmp_config["STATION_XML"] = tmp_xml
        update_stationXML(tmp_config)

if __name__ == "__main__":
    main()