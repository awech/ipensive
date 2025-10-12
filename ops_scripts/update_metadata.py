#!/home/seismicd/miniconda3/envs/ipensive/bin/python
import argparse
from obspy import UTCDateTime as utc
import ipensive.ipensive_utils as utils
from ipensive.metadata_utils import update_stationXML

if __name__ == "__main__":
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
    config, _ = utils.load_config(config_file)

    utils.setup_logging(utc.utcnow(), config)
    update_stationXML(config)