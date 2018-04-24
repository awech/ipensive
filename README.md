# web-infrasound


### Usage
edit <i>.config.py</i> with preprocessing parameters and array channel information<br>
and save as <i>config.py</i><br>
<br>
Run on a cron every 10 minutes:<br>
\*/10 \* \* \* \* python /path_to_file/array_processing.py >> /dev/null 2>&1<br>

### webpage configuration
Modify lines 545 - 561 (and 567) in <i>index.html</i> to match the networks and arrays in <i>config.py</i><br>
Make sure you have write privilages to 'out_dir' in <i>config.py</i> and that this directory is accessible from the web
<br><br>
(Note that <i>index.html</i> was originally written by Tom Parker)