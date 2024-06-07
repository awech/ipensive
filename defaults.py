"""
This file contains all default values and documentation for iPensive.
"""


def default_config():
    # Define default values for ipensive processing

    default_params = {

        ## SYSTEM PARAMETERS

        # CLIENT_TYPE Any Client understood by ObsPy
        # The default config uses an IP and port to configure an Earthworm Client. Other options include:
        #
        # 'CLIENT_TYPE': 'fdsn'
        # 'HOSTNAME': 'IRIS'
        # # 'PORT': < not used, can be omitted from conifg.py >
        #
        # An IRIS Client can also be defined this way for backward compatability:
        # # CLIENT_TYPE': < not used >
        # 'HOSTNAME': 'IRIS'
        # # 'PORT': < not used >
        #
        # A local FDSN web-service
        # 'CLIENT_TYPE': 'fdsn'
        # 'HOSTNAME': 'https://127.0.0.1'
        # 'PORT':18000
        #
        # A miniseed archive in SDS format
        # (https://docs.gempa.de/seiscomp/current/base/concepts/waveformarchives.html)
        # 'CLIENT_TYPE': 'sds'
        # 'HOSTNAME': '/path/to/miniseed/archive'  # path to the top level SDS directory
        # # 'PORT': < not used >
        #
        #
        'CLIENT_TYPE': "earthworm",
        'HOSTNAME': "pubavo.wr.usgs.gov",
        'PORT': 16022,

        # Logs and Output directories
        'OUT_WEB_DIR': '/www/avosouth.wr.usgs.gov/htdocs/infrasound',  # html & image ouptut directory
        'OUT_ASCII_DIR': '/www/avosouth.wr.usgs.gov/htdocs/infrasound/ascii_output', # ascii output directory (delete if undesired)
        'LOGS_DIR': '',

        # Extra links
        'EXTRA_LINKS': [
            {"AK Array Locations": "https://ds.iris.edu/gmap/#network=AV&location=01&channel=HDF&planet=earth"},
            {"AK Travel Times": "ak_traveltimes.html"},
            {"CNMI Array Locations": "https://ds.iris.edu/gmap/#network=MI,IM&station=FLX,H11N1,H11S1&planet=earth"},
            {"CNMI Travel Times": "cnmi_traveltimes.html"},
        ],  # added links for reference on the main ipensive page. Variable can be deleted if desired.

        # PROCESSING PARAMETERS
        'DURATION': 600,  # Duration of plotting image in seconds
        'LATENCY': 60,  # Seconds to wait before running array processing codes (allows all data to arrive)
        'EXTRA_PAUSE': 0,
        'FREQMIN': 0.8,  # lower filter band
        'FREQMAX': 5.0,  # upper filter band
        'TAPER': 5.0,  # taper length in seconds applied to data before processing
        'WINDOW_LENGTH': 30,  # window length in seconds for processing
        'OVERLAP': 15,  # overlap in seconds for processing windows
        'MIN_CHAN': 3,  # minimum number of channels to execute array processing codes

        # PLOTTING PARAMETERS
        'MCTHRESH': 0.6,  # minimum mcthreshold for plotting purposes
        'AZ_MIN': 0,  # minimum y-axis value for backazimuth panel
        'AZ_MAX': 360,  # maximum y-axis value for backazimuth panel
        'VEL_MIN': 0.25,  # minimum velocity for apaarent speed of sound box on velocity panel
        'VEL_MAX': 0.45,  # maximum velocity for apparent speed of sound box on velocity panel
        'ARRAY_LABEL': 'Infrasound',  # Label for the plots
    }

    return default_params
