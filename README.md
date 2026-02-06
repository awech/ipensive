# iPensive
Infrasound array processing for volcano monitoring. Developed for operational use at the Alaska Volcano Observatory and slightly generalized for broader use. The code processes infrasound and hydroacoustic array data in 10-minute increments, generating images and a webpage for displaying/navigating those images. Much of the package is a wrapper around Jordan Bishop's least trimmed squares [package](https://uaf-lts-array.readthedocs.io/en/master/index.html#), which performs the calculations to estimate trace velocity, back-azimuth, cross-correlation maxima, and flagged array element pairs (details: [https://doi.org/10.1093/gji/ggaa110](https://doi.org/10.1093/gji/ggaa110))

## Dependencies
![Static Badge](https://img.shields.io/badge/3.10%20%E2%80%93%203.13-blue?label=Python)

- obspy
- pandas
- jinja2
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
    conda create -n ipensive python=3.13
    conda activate ipensive
    ```

3. **Install ipensive as a package:**
    ```bash
    pip install .
    ```

    >Note: This will also install Jordan Bishop's [lts_array](https://github.com/uafgeotools/lts_array) which is included in ipensive package

## Usage 
#### Details and documentation found [here](/docs/index.md)

## Citation
This code may be cited directly as:
> Wech, A.W. (2026) ipensive (Version 1.0.1), U.S. Geological Survey Software Release, https://doi.org/10.5066/P14GHLCY.

The least trimmed squares infrasound methodology can be cited as:
>Bishop, J.W., Fee, D., & Szuberla, C. A. L., (2020). Improved infrasound array processing with robust estimators, Geophys. J. Int., 221(3) p. 2058-2074 doi: https://doi.org/10.1093/gji/ggaa110. 

And the accompanying `lts_array` code can be found here:
>https://github.com/uafgeotools/lts_array

## License and Disclaimer
[License](https://code.usgs.gov/vsc/seis/tools/ipensive/-/blob/main/LICENSE.md):
This project is in the public domain.

[Disclaimer](https://code.usgs.gov/vsc/seis/tools/ipensive/-/blob/main/DISCLAIMER.md):
This software is preliminary or provisional and is subject to revision.