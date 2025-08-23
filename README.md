# iPensive
Infrasound array processing for volcano monitoring. Developed for operational use at the Alaska Volcano Observatory and slightly generalized for broader use. The code processes infrasound and hydroacoustic array data in 10-minute increments, generating images and a webpage for displaying/navigating those images. Much of the package is a wrapper around Jordan Bishop's least trimmed squares [package](https://uaf-lts-array.readthedocs.io/en/master/index.html#), which performs the calculations to estimate trace velocity, back-azimuth, cross-correlation maxima, and flagged array element pairs (details: [https://doi.org/10.1093/gji/ggaa110](https://doi.org/10.1093/gji/ggaa110))

## Dependencies
![Static Badge](https://img.shields.io/badge/3.10%20%7C%203.11%20%7C%203.12-blue?label=Python)

- obspy
- pandas
- jinja2
- matplotlib
- pyyaml
- numba


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
    pip install -e .
    ```

    >Note: This will also install Jordan Bishop's [lts_array](https://github.com/uafgeotools/lts_array) which is included in ipensive package

## Usage 
#### Details and documentation found [here](/docs/index.md)