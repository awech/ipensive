"""
Regenerate test reference files after processing pipeline changes.

Produces:
  - tests/test_raw.mseed
  - tests/test_preprocessed.mseed
  - tests/test_results.csv

Requires network access (fetches data from FDSN/Earthscope).

Usage:
    python tests/generate_test_data.py
"""

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

from pathlib import Path
from obspy import UTCDateTime as utc
from ipensive import ipensive_utils as utils
from ipensive import data_utils, metadata_utils
from ipensive import array_processing as ap


IPENSIVE_CONFIG_PATH = Path("config/ipensive_config.yml")
T0 = "2025-08-19 20:30"
ARRAY = "Kenai"


def main():
    # Load config
    config = utils.load_ipensive_config(IPENSIVE_CONFIG_PATH)
    array_params = config[ARRAY]

    # Calculate time windows
    t1 = utc(T0) - array_params["DURATION"]
    t2 = utc(T0)
    T1 = t1 - array_params["TAPER"]
    T2 = t2 + array_params["WINDOW_LENGTH"] + array_params["TAPER"]

    # Grab raw data
    print("Grabbing raw data from FDSN...")
    st = data_utils.grab_data(array_params["CLIENT"], array_params["NSLC"], T1, T2)
    st.write("tests/test_raw.mseed", format="MSEED")
    print(f"  Wrote tests/test_raw.mseed ({len(st)} traces)")

    # Add metadata and preprocess
    st, *_ = metadata_utils.add_metadata(st, config, ARRAY, [])
    print("Preprocessing data...")
    st = data_utils.preprocess_data(st, t1, t2, array_params)
    st.write("tests/test_preprocessed.mseed", format="MSEED")
    print(f"  Wrote tests/test_preprocessed.mseed ({len(st)} traces)")

    # Run LTS to generate results CSV
    print("Running LTS analysis...")
    st, lat_list, lon_list = metadata_utils.add_metadata(st, config, ARRAY, [])
    array_params = utils.get_target_backazimuth(st, array_params)
    results_df, _ = ap.do_LTS(st, array_params, lat_list, lon_list, [])

    # Convert matplotlib date floats to ISO time strings for the CSV
    from matplotlib import dates
    results_df["Time"] = [dates.num2date(ti).strftime('%Y-%m-%dT%H:%M:%S') for ti in results_df["Time"]]
    results_df.to_csv("tests/test_results.csv", index=False, float_format="%.7f")
    print(f"  Wrote tests/test_results.csv ({len(results_df)} rows)")

    print("\nDone.")


if __name__ == "__main__":
    main()
