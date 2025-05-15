# iPensive
Infrasound array processing for volcano monitoring. Developed for operational use at the Alaska Volcano Observatory and slightly generalized for broader use.

### Dependencies
Tested on Python 3.12

- obspy
- pandas
- jinja2
- matplotlib
- pyyaml
- lts_array (from GitHub: [uafgeotools/lts_array](https://github.com/uafgeotools/lts_array))

### Installation

1. **Clone this repository:**
    ```bash
    git clone https://github.com/awech/ipensive.git
    cd ipensive
    ```

2. **(Recommended) Create and activate a virtual environment:**
    ```bash
    conda create -n ipensive python=3.12
    conda activate ipensive
    ```

3. **Install all dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

    This will also install [lts_array](https://github.com/uafgeotools/lts_array) directly from GitHub.

### Configuration
Edit [data_source.yml](data_source.yml) with the obspy client type (e.g., *earthworm, fdsn, sds, local_fdsn*) and relevant parameters (e.g., *HOSTNAME, PORT, DIRECTORY, etc.*).

Edit [config.yml](config.yml). 
1. **system  parameters**: defines input files
    - DATA_SOURCE: [data_source.yml](data_source.yml) - data source parameters
    - STATION_XML: [station.xml](station.xml) - station metadata
    - TARGETS_FILE: [volcano_list.txt](volcano_list.txt) - csv file of locations with backazimuths of interest (Target,Longitude,Latitude)
    - OUT_WEB_DIR - location to write out html and images
    - OUT_ASCII_DIR - location to write results output in ascii format
    - LOGS_DIR - location of log output when run as a cron
2. **processing parameters**
    - processing parameters: controlling data processing details (filters, window length, etc.)
    - plotting parameters: slight control on how a few things appear (mostly this allows for differentiating between acoustic and hydroacoustic velocities)
    - note, these are defaults for all arrays, but each individual array can have its own unique parameters to selectively override the default
3. **network parameters**
    - defines network and array structure for the web output
4. **array parameters**
    - NSLC: list of array channels. Metadata must be in ```STATION_XML```, or can include lat/lon/gain information manually here
    - TARGETS: list of targets for which you want to plot back-azimuths. Must be in ```TARGETS_FILE```, or can include lat/lon information manually here 


### Usage
1. **Automatically on a cron**

    ```*/10 * * * * /path_to_file/array_processing.py >> /dev/null 2>&1```
2. **Manually**

    ```python array_processing.py -c <config_file> -t <yyyymmddHHMM>```

    see ```python array_processing.py -h``` for more options. In particular, the ```-a``` option is useful for processing a single array if more than one are defined in the config file.


### Webpage configuration
The webpage is automatically generated each time using jinja2 to populate the [index.template](index.template) file with the network and array structure from the config file.
(This method was adapted from a forked version of this code by Israel Brewster, and the html and javascript front-end was originally developed by Tom Parker.)

### Back populate data
Modify T1 and T2 in [back_populate.py](back_populate.py) and run this script to generate images between time T1 and T2. There are options to overwrite existing files (or not) or only process a subset of arrays.