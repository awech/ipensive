# iPensive
Infrasound array processing for volcano monitoring. Developed for operational use at the Alaska Volcano Observatory and slightly generalized for broader use.

### Dependencies
Tested on python 3.7<br>

- obspy<br>
- pandas<br>
- jinja2<br>

### Usage
Edit <i>config.py</i>. There are<br>
- system level parameters: controlling data source and output location<br>
- processing parameters: controlling data processing details (filters, window length, etc.)<br>
- plotting parameters: slight control on how a few things appear (mostly this allows for differentiating between acoustic and hydroacoustic velocities)<br>

Note that each individual array can have its own unique parameters to selectively override the default (though this wouldn't work well for the output directories)<br>
<br>
Run on a cron every 10 minutes:<br>
`\*/10 \* \* \* \* python /path_to_file/array_processing.py >> /dev/null 2>&1` <br>


### webpage configuration
The webpage is automatically generated each time using jinja2 to populate the <i>index.template</i> file with the network and array structure from the config file.
(This method was adapted from a forked version of this code by Israel Brewster, and the html and javascript front-end was originally developed by Tom Parker)<br>
Note that if you do not want to process a particular array (e.g. bad/no data), the array must remain defined in <i>config.py</i> in order to remain a viewing option on the webpage (for past data). In this case, you could just comment out the SCNL's for that particular array in <i>config.py</i> 

### Back populate data
Modify T1 and T2 in <i>back_populate.py</i> and run this script to generate images between time T1 and T2