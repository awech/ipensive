## [1.0.0] - 2025-08-26

**Not Yet Approved Release**

This version is the official public release of ipensive and the start of semantic versioning. This version is a massive overhaul from the original Github code base, which is included as a release here as version 0.0.1. Changes include how ipensive is installed, configured and run. The full commit history is preserved in this repository. Notable changes are listed below:

### Changed
- moved scripts to `ops_scripts`
- added command line arguments to scripts
- structured as package for simpler installation
- uses Jordan Bishop's [lts_array](https://github.com/uafgeotools/lts_array) package
- move LTS code and ipensive code to `src/` subdirectory
- improved documentation
- modified for python >= 3.10
- implement `matplotlib`'s `subplot_mosaic`
- change `.py` config files to  `.yml` files
- moved channel metadata out of array config file and into `.xml` file
- moved infrasound target locations to separate `.csv` file
- separate array configuration from system configuration
- added support for various `obspy` Clients
- implement pytests
- added Dockerfile