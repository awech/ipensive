iPensive is a tool for automatic processing of infrasound array data for
volcano monitoring. It creates .png images as results and stores them in a
way that can be visualized through an html file, index.html. Array
processing results can also be written to an ASCII file for additional
processing.

## Installation

Download the [latest release](https://github.com/awech/ipensive/tree/master)
or use `git` to clone the entire repository to a working directory.

iPensive was developed with Python 3.7 and uses the following dependencies:  
- [obspy](http://www.obspy.org/)
- [pandas](http://pandas.pydata.org/)
- [jinja2](https://palletsprojects.com/p/jinja/) 

These dependencies can be installed via [Anaconda](https://www.anaconda.com/download)
on the command line. An environment.yml file does not exist at this time, but
creating the environment and installing the dependencies can be done on the command
line.

Assuming you have installed Anaconda, create the ipensive environment with
Python 3.7
```
conda env create -n ipensive python=3.7
```

```
conda config --add channels conda-forge
conda install -c conda-forge obspy pandas jinja2
```

Finally, activate the environment:
```
conda activate ipensive
```

Later, you will need to know the path to the Python executable for this
environment. You can find it with the following command:
```
which python
```

## Configuration

Edit <i>config.py</i>. There are
- system level parameters: controlling data source and output location
- processing parameters: controlling data processing details (filters, window length, etc.)
- plotting parameters: slight control on how a few things appear (mostly this allows for differentiating between acoustic and hydroacoustic velocities)

Default values for all config variables are defined in a dictionary in
ipensive_utils.default_config(). An explanation of each variable is also
provided in this file.

Note that each individual array can have its own unique parameters to
selectively override the default (though this wouldn't work well for the
output directories)

## Usage

iPensive can be run on the command line, or it can be set up as cronjob.
This example shows how to set up iPensive as a cron every 10 minutes:
```
*/10 * * * * /path_to_file/array_processing.py >> /dev/null 2>&1
```