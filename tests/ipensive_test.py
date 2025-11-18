# pip install pytest
# pip install coverage
# To test: `pytest -v --disable-warnings`
# -v for verbose
# --disable-warnings (does what it says)
# To check coverage:
# `coverage run -m pytest` to collect coverage data
# `coverage html -d ~/Desktop/test_html` to generate HTML report

import os

from ipensive import ipensive_utils as utils
from ipensive import data_utils, metadata_utils, plotting_utils
from ipensive import array_processing as ap

CURR_DIR = os.path.dirname(__file__)
IPENSIVE_CONFIG_PATH = os.path.abspath(os.path.join(CURR_DIR, '../config/example_1/ipensive_config.yml'))
ARRAYS_CONFIG_PATH = os.path.abspath(os.path.join(CURR_DIR, '../config/example_1/arrays_config.yml'))
os.chdir(CURR_DIR)

T0 = "2025-08-19 20:30"
ARRAY = "Kenai"


def test_parse_args_defaults(monkeypatch):
    # Patch sys.argv
    import sys
    monkeypatch.setattr(sys, 'argv', ['array_processing.py'])
    args = ap.parse_args()
    assert hasattr(args, 'config')


def test_starttime(monkeypatch):
    import sys
    from obspy import UTCDateTime

    T1 = UTCDateTime(T0).strftime("%Y%m%d%H%M")

    monkeypatch.setattr(sys, 'argv', ['array_processing.py', '--time', T1])
    args = ap.parse_args()

    config = utils.load_ipensive_config(IPENSIVE_CONFIG_PATH)    
    T1, delay = ap.get_starttime(config, args)

    assert T0 == UTCDateTime(T1)


def test_ipensive_logger_output(monkeypatch, tmp_path, caplog):
    import sys
    import logging

    monkeypatch.setattr(sys, 'argv', ['array_processing.py'])
    args = ap.parse_args()

    monkeypatch.setenv("FROMCRON", "yep")
    
    config = utils.load_ipensive_config(IPENSIVE_CONFIG_PATH)
    config["LOGS_DIR"] = tmp_path
    T0, delay = ap.get_starttime(config, args)

    utils.setup_logging(T0, config)
    logger = logging.getLogger("ipensive")
    test_message = "Logger output test message"
    with caplog.at_level(logging.INFO, logger="ipensive"):
        logger.info(test_message)
    assert any(test_message in m for m in caplog.messages)

    assert os.path.exists(tmp_path / T0.strftime("%Y") / T0.strftime("%Y-%m"))


def test_load_config_and_arrays():
    config = utils.load_ipensive_config(IPENSIVE_CONFIG_PATH)
    assert 'NETWORKS' in config
    assert 'array_list' in config
    assert isinstance(config['array_list'], list)
    assert 'OUT_WEB_DIR' in config
    assert 'OUT_ASCII_DIR' in config
    assert 'LOGS_DIR' in config
    # Check that arrays from arrays_config are present
    for net in config['NETWORKS']:
        for arr in config['NETWORKS'][net]:
            assert arr in config


def test_env_variable_load_config(monkeypatch, tmp_path):
    import shutil

    env_ipensive_file = str(tmp_path / "ipensive_config_env.yml")
    env_array_file = str(tmp_path / "arrays_config_env.yml")
    shutil.copy(IPENSIVE_CONFIG_PATH, env_ipensive_file)
    shutil.copy(ARRAYS_CONFIG_PATH, env_array_file)

    monkeypatch.setenv("IPENSIVE_CONFIG", env_ipensive_file)
    monkeypatch.setenv("ARRAYS_CONFIG", env_array_file)

    ipensive_config_file = utils.get_config_file()
    config2 = utils.load_ipensive_config(ipensive_config_file)

    assert str(ipensive_config_file) == env_ipensive_file
    assert utils.check_path(config2["array_config_files"][0]) == utils.check_path(env_array_file)


def test_get_pngfile_path():
    from pandas import Timestamp
    config = utils.load_ipensive_config(IPENSIVE_CONFIG_PATH)
    t = Timestamp("2024-01-01T12:00:00")
    arr = config['array_list'][0]
    path = utils.get_pngfile_path(t, arr, config)
    assert str(path).endswith('.png')


def test_get_obspy_client_fdsn():
    # Minimal config for FDSN
    conf = {'CLIENT_TYPE': 'fdsn', 'HOSTNAME': 'IRIS', 'TIMEOUT': 5}
    client = utils.get_obspy_client(conf)
    assert hasattr(client, 'get_waveforms')


def test_update_metadata(tmp_path):
    from obspy import read
    
    config = utils.load_ipensive_config(IPENSIVE_CONFIG_PATH)
    config["STATION_XML"] = tmp_path / "station.xml"

    st = read(f"{CURR_DIR}/test_raw.mseed")
    assert metadata_utils.check_FDSN(st[0], config[ARRAY]["CLIENT"])

    metadata_utils.update_stationXML(config)
    assert config["STATION_XML"].exists()


def test_data_and_preprocessing():

    from obspy import read
    from obspy import UTCDateTime as utc
    from numpy.testing import assert_allclose


    config_file = utils.get_config_file()
    config = utils.load_ipensive_config(config_file)
    array_params = config[ARRAY]

    t1 = utc(T0) - array_params["DURATION"]
    t2 = utc(T0)

    T1 = t1 - array_params["TAPER"]  # Extended start time for tapering
    T2 = t2 + array_params["WINDOW_LENGTH"] + array_params["TAPER"]  # Exte

    st = data_utils.grab_data(array_params["CLIENT"], array_params["NSLC"], T1, T2)
    ST = read(f"{CURR_DIR}/test_raw.mseed")

    for tr in st:
        TR = ST.select(id=tr.id)[0]
        assert_allclose(tr.data, TR.data, atol=1e-5, rtol=1e-8)

    st, *_ = metadata_utils.add_metadata(st, config, ARRAY, [])
    st = data_utils.preprocess_data(st, t1, t2, [], array_params)
    ST = read(f"{CURR_DIR}/test_preprocessed.mseed")

    for tr in st:
        TR = ST.select(id=tr.id)[0]
        assert_allclose(tr.data, TR.data, atol=1e-9, rtol=1e-12)


def test_write_html(tmp_path):

    config = utils.load_ipensive_config(IPENSIVE_CONFIG_PATH)
    # Set a temporary output directory for HTML
    config['OUT_WEB_DIR'] = tmp_path
    # Add required keys for template rendering if missing
    if 'EXTRA_LINKS' not in config:
        config['EXTRA_LINKS'] = []

    # Run write_html
    utils.write_html(config)
    # Check that the HTML file was created
    html_file = config['OUT_WEB_DIR'] / 'index.html'
    assert html_file.exists()
    with open(html_file) as f:
        content = f.read()
    assert '<html' in content or 'DOCTYPE html' in content or 'body' in content


def test_write_data_files(tmp_path):
    import pandas as pd
    from obspy import UTCDateTime, read
    
    config = utils.load_ipensive_config(IPENSIVE_CONFIG_PATH)
    config['OUT_ASCII_DIR'] = tmp_path / "ascii"
    config['OUT_VALVE_DIR'] = tmp_path / "valve"

    DF = pd.read_csv(f"{CURR_DIR}/test_results.csv")
    st = read(f"{CURR_DIR}/test_preprocessed.mseed")
    utils.write_data_files(UTCDateTime(T0), st, DF, config)

    # Check files exist
    found = list(tmp_path.rglob('*.txt'))
    assert len(found) == 2
    for f in found:
        df = pd.read_csv(f, sep=None, engine='python')
        assert not df.empty


def test_LTS_and_image_output(tmp_path):
    from obspy import read
    from obspy import UTCDateTime as utc
    from pandas import read_csv
    from pandas.testing import assert_frame_equal

    config = utils.load_ipensive_config(IPENSIVE_CONFIG_PATH)
    config['OUT_WEB_DIR'] = tmp_path

    t2 = utc(T0)
    array_params = config[ARRAY]
    skip_chans = []
    st = read(f"{CURR_DIR}/test_preprocessed.mseed")

    st, lat_list, lon_list = metadata_utils.add_metadata(st, config, ARRAY, skip_chans)
    array_params = utils.get_target_backazimuth(st, array_params)
    test_df, lts_dict = ap.do_LTS(st, array_params, lat_list, lon_list, skip_chans)

    expected_df = read_csv(f"{CURR_DIR}/test_results.csv")

    assert_frame_equal(expected_df, test_df.round(7), rtol=1.5e-5)

    utils.web_folders(t2, config, array_params)
    for plotsize in ["big", "thumbnail"]:
        fig = plotting_utils.plot_results(
            t2, st, test_df, lts_dict, skip_chans, array_params, plotsize
        )
        plotting_utils.save_figure(fig, config, array_params, t2, plotsize)

    # Check that the big image file was created
    image_files = list(tmp_path.rglob('*0.png'))
    assert len(image_files) > 0
    for f in image_files:
        assert f.exists()

    # Check that the thumbnail image created
    image_files = list(tmp_path.rglob('*_thumb.png'))
    assert len(image_files) > 0
    for f in image_files:
        assert f.exists()


def test_array_processing(tmp_path):
    from obspy import UTCDateTime as utc
    config = utils.load_ipensive_config(IPENSIVE_CONFIG_PATH)

    config['OUT_ASCII_DIR'] = tmp_path / "ascii"
    config['OUT_VALVE_DIR'] = tmp_path / "valve"
    config['OUT_WEB_DIR'] = tmp_path / "web"

    config["plot"] = True

    fig = None
    fig = ap.process_array(config, ARRAY, utc(T0), return_figure=True)

    assert fig is not None

    image_files = list(tmp_path.rglob('*.png'))
    assert len(image_files) == 0


def test_manual_metadata(tmp_path):
    from obspy import UTCDateTime as utc
    config = utils.load_ipensive_config(IPENSIVE_CONFIG_PATH)

    config['OUT_ASCII_DIR'] = tmp_path / "ascii"
    config['OUT_VALVE_DIR'] = tmp_path / "valve"
    config['OUT_WEB_DIR'] = tmp_path / "web"

    config["plot"] = True

    ap.process_array(config, "Wake Island South", utc(T0), return_figure=False)

    image_files = list(tmp_path.rglob('*.png'))
    assert len(image_files) == 2


def test_multiple_array_configs(tmp_path):
    
    config1 = utils.load_ipensive_config(IPENSIVE_CONFIG_PATH)
    IPENSIVE_CONFIG_PATH2 = os.path.abspath(os.path.join(CURR_DIR, '../config/example_2/ipensive_config.yml'))
    config2 = utils.load_ipensive_config(IPENSIVE_CONFIG_PATH2)

    assert config1.keys() == config2.keys()

    for arr in config1['array_list']:
        assert arr in config2['array_list']
        assert config1[arr]["STATION_XML"] != config2[arr]["STATION_XML"]
        assert config1[arr]["TARGETS_FILE"] != config2[arr]["TARGETS_FILE"]

    assert len(config1["array_config_files"]) == 1
    assert len(config2["array_config_files"]) == 2