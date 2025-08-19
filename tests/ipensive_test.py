import os
import types

# Third-party and package modules
import ipensive.ipensive_utils as utils
import ipensive.array_processing as ap

CONFIG_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '../config/ipensive_config.yml'))
ARRAYS_CONFIG_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '../config/arrays_config.yml'))

## TODO
# Add tests for the new features and improvements
# Cleant this up

# TODO
# test 1) pre-processing and 2) array processing
# with raw, cut seed data
def test_preprocessing_data():
    from obspy import UTCDateTime as utc
    T1 = utc("2024-01-01T12:00:00")
    T2 = utc("2024-01-01T12:00:00")

# TODO
# Test image output

def test_get_config_file():
    # Should return the default config path
    config_file = utils.get_config_file()
    assert str(config_file).endswith('ipensive_config.yml')


def test_load_config_and_arrays():
    config = utils.load_config(CONFIG_PATH)
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


def test_env_variable_load_config(monkeypatch):
    config1 = utils.load_config(CONFIG_PATH)
    monkeypatch.setenv("IPENSIVE_CONFIG", str(CONFIG_PATH))
    monkeypatch.setenv("ARRAYS_CONFIG", str(ARRAYS_CONFIG_PATH))
    config_file = utils.get_config_file()
    assert config1 == utils.load_config(config_file)


def test_get_obspy_client_fdsn():
    # Minimal config for FDSN
    conf = {'CLIENT_TYPE': 'fdsn', 'HOSTNAME': 'IRIS', 'TIMEOUT': 5}
    client = utils.get_obspy_client(conf)
    assert hasattr(client, 'get_waveforms')


def test_get_file_path():
    from pandas import Timestamp
    config = utils.load_config(CONFIG_PATH)
    t = Timestamp("2024-01-01T12:00:00")
    arr = config['array_list'][0]
    path = utils.get_file_path(t, arr, config)
    assert str(path).endswith('.png')


def test_write_html(tmp_path):
    from pathlib import Path
    config = utils.load_config(CONFIG_PATH)
    # Set a temporary output directory for HTML
    config['OUT_WEB_DIR'] = str(tmp_path)
    # Add required keys for template rendering if missing
    if 'network_list' not in config:
        config['network_list'] = list(config['NETWORKS'].keys())
    if 'EXTRA_LINKS' not in config:
        config['EXTRA_LINKS'] = []

    # Run write_html
    ap.write_html(config)
    # Check that the HTML file was created
    html_file = Path(config['OUT_WEB_DIR']) / 'index.html'
    assert html_file.exists()
    with open(html_file) as f:
        content = f.read()
    assert '<html' in content or 'DOCTYPE html' in content or 'body' in content


def test_write_ascii_file_and_valve_file(tmp_path):
    import numpy as np
    import pandas as pd
    from obspy import UTCDateTime
    config = utils.load_config(CONFIG_PATH)
    arr = config['array_list'][0]
    t2 = UTCDateTime("2024-01-01T12:00:00")
    t = ["2024-01-01 11:59:00", "2024-01-01 12:00:00"]
    pressure = np.array([1.0, 2.0])
    azimuth = np.array([10.0, 20.0])
    velocity = np.array([0.3, 0.4])
    mccm = np.array([0.8, 0.9])
    rms = np.array([0.1, 0.2])
    config['OUT_ASCII_DIR'] = str(tmp_path)
    config['OUT_VALVE_DIR'] = str(tmp_path)
    utils.write_ascii_file(t2, t, pressure, azimuth, velocity, mccm, rms, arr, config)
    utils.write_valve_file(t2, t, pressure, azimuth, velocity, mccm, rms, arr, config)
    # Check files exist
    found = list(tmp_path.rglob('*.txt'))
    assert found
    for f in found:
        df = pd.read_csv(f, sep=None, engine='python')
        assert not df.empty


def test_process_array_minimum(monkeypatch):
    # Test process_array with minimal config and mocks
    import numpy as np
    from obspy import UTCDateTime, Stream, Trace
    config = utils.load_config(CONFIG_PATH)
    arr = config['array_list'][0]
    # Patch utils.grab_data to return a valid Stream
    def fake_grab_data(client, NSLC, T1, T2, fill_value=0):
        tr = Trace()
        tr.id = NSLC[0] if isinstance(NSLC, list) else list(NSLC.keys())[0]
        tr.stats.sampling_rate = 100
        tr.stats.starttime = T1
        tr.data = np.ones(int((T2-T1)*100))
        st = Stream([tr])
        tr.stats.coordinates = types.SimpleNamespace(latitude=60.0, longitude=-150.0)
        tr.inventory = None
        return st
    monkeypatch.setattr(utils, 'grab_data', fake_grab_data)
    # Patch plotting and file writing
    monkeypatch.setattr(ap, 'plot_results', lambda *a, **k: None)
    monkeypatch.setattr(utils, 'web_folders', lambda *a, **k: None)
    monkeypatch.setattr(utils, 'write_valve_file', lambda *a, **k: None)
    monkeypatch.setattr(utils, 'write_ascii_file', lambda *a, **k: None)
    # Patch ltsva
    monkeypatch.setattr(ap, 'ltsva', lambda *a, **k: (np.array([0.3]), np.array([10.0]), np.array([1.0]), np.array([0.8]), {}, np.array([0.1]), None))
    # Patch get_target_backazimuth
    monkeypatch.setattr(utils, 'get_target_backazimuth', lambda st, config, arr_params: arr_params)
    # Patch add_metadata
    monkeypatch.setattr(utils, 'add_metadata', lambda st, config, skip_chans: st)
    # Patch add_coordinate_info
    monkeypatch.setattr(utils, 'add_coordinate_info', lambda st, config, arr: st)
    # Patch Stream.slide
    monkeypatch.setattr(Stream, 'slide', lambda self, window_length, step: [self])
    # Patch config['plot']
    config['plot'] = False
    T0 = UTCDateTime("2024-01-01T12:00:00")
    ap.process_array(config, arr, T0)


def test_parse_args_defaults(monkeypatch):
    # Patch sys.argv
    import sys
    monkeypatch.setattr(sys, 'argv', ['array_processing.py'])
    args = ap.parse_args()
    assert hasattr(args, 'config')
