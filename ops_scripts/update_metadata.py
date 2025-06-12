import argparse
from pathlib import Path
from obspy import UTCDateTime as utc
from ipensive.ipensive_utils import setup_logging, load_config
from ipensive.metadata import update_stationXML

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        epilog="e.g.: python update_metadata.py -c <filename.yml>"
    )

    default_config = Path(__file__).parent.parent / "config/config.yml"

    parser.add_argument(
        "-c",
        "--config",
        type=str,
        help="Name of the config file (yml)",
        default=default_config,
    )
    args = parser.parse_args()
    config_file = args.config
    config = load_config(config_file)

    setup_logging(utc.utcnow(), config)
    update_stationXML(config)