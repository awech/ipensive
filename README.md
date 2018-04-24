# web-infrasound


### Usage
edit <i>.config.py</i> with preprocessing parameters and array channel information<br>
and save as <i>config.py</i><br>
<br>
Run on a cron every 10 minutes:<br>
\*/10 \* \* \* \* python /path_to_file/array_processing.py >> /dev/null 2>&1<br>

### webpage configuration
Modify lines 545 - 561 in <i>index.html</i> to match the networks and arrays in <i>config.py</i><br>
Make sure that <i>index.html</i> is in directory as 'out_dir' in <i>config.py</i>