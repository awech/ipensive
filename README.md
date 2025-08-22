# iPensive
Infrasound array processing for volcano monitoring. Developed for operational use at the Alaska Volcano Observatory and slightly generalized for broader use.

## Dependencies
![Static Badge](https://img.shields.io/badge/3.10%20%7C%203.11%20%7C%203.12-blue?label=Python)


- obspy
- pandas
- jinja2
- matplotlib
- pyyaml
- lts_array (from GitHub: [uafgeotools/lts_array](https://github.com/uafgeotools/lts_array))

## Installation

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

## Quick-start
Change in the `/ops_scripts` subdirectory 
```bash
ipensive
├── config
├── data
├── ipensive
├── ops_scripts
│   ├── run_processing.py
│   ├── update_metadata.py
│   └── back_populate.py
├── templates
├── tests
└── ...
```

and execute the `run_processing.py` script:
```bash
cd ops_scripts
python run_processing.py
```

This should create an `/output` subdirectory with an `html` folder and `ascii_output` folder:
```bash
ipensive
├── config
├── ipensive
├── ops_scripts
...
├── output
│   ├── html
│   └── ascii_output
└── ...
```

Navigate to and open `/output/html/index.html` with your web browser to verify the webpage and images were generated.


## Configuration Files
There are 2 configuration files in `/config`:
```bash
ipensive
├── config
│   ├── ipensive_config.yml
│   └── arrays_config.yml
├── ipensive
├── ops_scripts
└── ...
```

1.  [`ipensive_config.yml`](../config/ipensive_config.yml) controls the data source and the output locations
2.  [`arrays_config.yml`](../config/arrays_config.yml) defines the arrays, processing parameters, and plot controls

Both files have example/template entries in them to demonstrate the structure and available configuration options. The simplest step would be to edit these files in place with your arrays, data source, and desired output locations. The path to a separate `ipensive_config.yml` can also be provided as an argument to the scripts in `/ops_scripts`, and the path to `arrays_config.yml` can also be defined within `ipensive_config.yml`.

>NOTE:
>There is also the option to create an environment variables `IPENSIVE_CONFIG` and `ARRAYS_CONFIG` defining paths to separate locations of these configs. 


### System config: `ipensive_config.yml`
This config file defines
1. Data Source:
    - `CLIENT_TYPE` (e.g., *earthworm, fdsn, sds, local_fdsn*)
    - `HOSTNAME`
    - `PORT` (if relevant)
    - `TIMEOUT` (optional. Defaults to 30s)

2. Output Directories:
    - `OUT_WEB_DIR`: <directory_path> location to output html and images
    - `OUT_ASCII_DIR`: <directory_path> location to output ascii files of processing results
    - `LOGS_DIR`: <directory_path> location to output log files (*Note: must set environment variable* `FROMCRON=yep` *to write log output*)

3. Metadata and back-azimuth target information:
    - `STATION_XML`: <path_to_file> (station.xml file updated routinely by `ops_sripts/update_metadata.py`)
    - `TARGETS_FILE`: <path_to_file> (.csv file with columns: `Target`,`Longitude`,`Latitude`)


### Array config: `arrays_config.yml`
This config file defines:

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
Operational scripts are located in `/ops_scripts` directory

1. **Manually**

    ```python run_processing.py.py -c <config_file> -t <yyyymmddHHMM>```

     `-c` and `-t` are optional. If `-c` is not supplied, the config defaults to either the path defined by environment varable `IPENSIVE_CONFIG` (1st) or `/config/ipensive_config.yml` (2nd). If `-t` is not supplied, the time is rounded to the most recent 10-minute mark. See ```python run_processing.py -h``` for more options. In particular, the ```-a``` option is useful for processing a single array if more than one are defined in the arrays config file.

2. **Automatically on a cron**

    ```*/10 * * * * python /path_to_file/run_processing.py.py >> /dev/null 2>&1```

3. **Back populate**

    ```bash
    python back_populate.py -s 202507010000 -e 202507020000
    ```

    see ```python back_populate.py -h``` for more options.

4. **Update Metadata**
    ```bash
    python update_metadata.py -c <config_file>
    ```
    This will update the `data/stations.xml` file with station-channel lat/lon information used for array processing. Ideally run somewhat periodically, though not as often as the array processing itself.


### Webpage configuration
The webpage is automatically generated each time using jinja2 to populate the [index.template](index.template) file with the network and array structure from the config file.
(This method was adapted from a forked version of this code by Israel Brewster, and the html and javascript front-end was originally developed by Tom Parker.)