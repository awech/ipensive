# Release Notes

## Versions

To check which version of `ipensive` you have installed:
```
pip show ipensive
```

## Maintenance team

Current and past maintainers of `ipensive`:
- [@awech](https://github.com/awech)

## [1.0.1] - 2026-02-06

**Minor Patch** - [https://doi.org/10.5066/P14GHLCY](https://doi.org/10.5066/P14GHLCY)

This release fixes a relative paths bug introduced when converting scripts to command-line executables. 

### :bug: Bugs
- fixed relative pathing issue in `ipensive_utils` when no default config file is set. This is mainly an issue that affects its usability right out of the box

## [1.0.0] - 2025-11-25

**Fist Approved Release** - [https://doi.org/10.5066/P1JBRDDF](https://doi.org/10.5066/P1JBRDDF)

This version is the official public release of ipensive and the start of semantic versioning. This version is a massive overhaul from the original Github code base, which is included as a release here as version 0.0.1. Changes include how ipensive is installed, configured and run. The full commit history is preserved in this repository. Notable changes are listed below:


### :sparkles: Features
- uses Jordan Bishop's [lts_array](https://github.com/uafgeotools/lts_array) package
- added support for various `obspy` Clients
- modified for python >= 3.10. Tested on 3.10-3.13
- support for multiple sets of array configurations and targets
- added command line arguments to scripts

### :books: Documentation
- improved documentation
- implement pytests
- improved logging with timestamps

### :gear: Configuration
- change `.py` config files to  `.yml` files
- moved channel metadata out of array config file and into `.xml` file
- moved infrasound target locations to separate `.csv` file
- separated array configuration from system configuration
- new option to use environment variables pointing to config files

### :hammer_and_wrench: Under the hood
- structured as package for simpler installation
- implement `matplotlib`'s `subplot_mosaic`
- added Dockerfile for containerization
- move LTS code and ipensive code to `src/` subdirectory

### :bug: Bugs
- probably