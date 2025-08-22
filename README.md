# iPensive
Infrasound array processing for volcano monitoring. Developed for operational use at the Alaska Volcano Observatory and slightly generalized for broader use.

### Dependencies
![Static Badge](https://img.shields.io/badge/3.10%20%7C%203.11%20%7C%203.12-blue?label=Python)


- obspy
- pandas
- jinja2
- matplotlib
- pyyaml
- lts_array (from GitHub: [uafgeotools/lts_array](https://github.com/uafgeotools/lts_array))

### Installation

1. **Clone this repository:**
    ```bash
    git clone https://code.usgs.gov/vsc/seis/tools/ipensive.git
    cd ipensive
    ```

2. **(Recommended) Create and activate a virtual environment:**
    ```bash
    conda create -n ipensive python=3.12
    conda activate ipensive
    ```

3. **Install ipensive as a package:**
    ```bash
    pip install -r requirements.txt
    pip install -e .
    ```

    This will also install [lts_array](https://github.com/uafgeotools/lts_array) directly from GitHub.

### System Configuration
Edit [ipensive_config.yml](../config/ipensive_config.yml) in `../config` or create an environment variable `IPENSIVE_CONFIG` defining the location of the main config yml file. This config file defines the Data Source:
- `CLIENT_TYPE` (e.g., *earthworm, fdsn, sds, local_fdsn*)
- `HOSTNAME`
- `PORT` (if relevant)
- `TIMEOUT` (optional. Defaults to 30s)

Output Directories:
- `OUT_WEB_DIR`: <directory_path> location to output html and images
- `OUT_ASCII_DIR`: <directory_path> location to output ascii files of processing results
- `LOGS_DIR`: <directory_path> location to output log files (*Note: must set environment variable* `FROMCRON=yep` *to write log output*)

Metadata and back-azimuth target information:
- `STATION_XML`: <path_to_file> (station.xml file updated routinely by `ops_sripts/update_metadata.py`)
- `TARGETS_FILE`: <path_to_file> (.csv file with columns: `Target`,`Longitude`,`Latitude`)


### Array Configuration
Edit [arrays_config.yml](../config/arrays_config.yml) in `../config` or create an environment variable `ARRAYS_CONFIG` defining the location of the arrays config yml file. This config file defines:

1. **processing parameters** `PARAMS`
    - processing parameters: controlling data processing details (filters, window length, etc.)
    - plotting parameters: slight control on how a few things appear (mostly this allows for differentiating between acoustic and hydroacoustic velocities)
    - *NOTE*: these are defaults for all arrays, but each individual array can have its own unique parameters to selectively override the default
3. **network parameters** `NETWORKS`
    - network and array structure for the web output
4. **array parameters** `<ARRAY_NAME>`
    - `NSLC`: list of array channels. Metadata must be in `STATION_XML` file, or can include lat/lon/gain information manually here
    - `TARGETS`: list of targets for which you want to plot back-azimuths. Must be in `TARGETS_FILE`, or can include lat/lon information manually here
    - Any parameter set in `PARAMS` above can be customized for individual arrays here


### Usage
Operational scripts are located in `../ops_scripts` directory

1. **Automatically on a cron**

    ```*/10 * * * * python /path_to_file/run_processing.py.py >> /dev/null 2>&1```

2. **Manually**

    ```python run_processing.py.py -c <config_file> -t <yyyymmddHHMM>```

    see ```python run_processing.py -h``` for more options. In particular, the ```-a``` option is useful for processing a single array if more than one are defined in the config file.

3. **Back populate**
    - ```python back_populate.py -s 202507010000 -e 202507020000```

    see ```python back_populate.py -h``` for more options.


### Webpage configuration
The webpage is automatically generated each time using jinja2 to populate the [index.template](index.template) file with the network and array structure from the config file.
(This method was adapted from a forked version of this code by Israel Brewster, and the html and javascript front-end was originally developed by Tom Parker.)