#!/Users/awech/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
Created on 24-Apr-2018
@author: awech
"""

import os
import sys
import numpy as np
from pandas import DataFrame
from obspy.core import UTCDateTime
from matplotlib import dates
import time
import config
import utils
import warnings
os.chdir(config.working_dir)


if len(sys.argv) == 1:                              # no time given, use current time
    T0=UTCDateTime.utcnow()                         # get current timestamp
    T0=UTCDateTime(T0.strftime('%Y-%m-%d %H:%M')[:-1]+'0') # round down to the nearest 10-minute
    time.sleep(config.latency+config.window_length)
else:                                               # time given, use it
    if len(sys.argv)==2:                            # time given as single string (eg. 201705130301)
        T0 = sys.argv[1]
    elif len(sys.argv)==3:                          # time given as 2 strings (eg. 20170513 03:01)
        T0='{}{}'.format(sys.argv[1],sys.argv[2])
    else:
        warnings.warn('Too many input arguments. eg: array_processing.py 201701020205')
        sys.exit()      
    try:
        T0 = UTCDateTime(T0)
        T0 = UTCDateTime(T0.strftime('%Y-%m-%d %H:%M')[:-1]+'0') # round down to the nearest 10-minute
    except:
        warnings.warn('Needs end-time argument. eg: array_processing.py 201701020205')
        sys.exit()


t1 = T0-config.duration
t2 = T0
print('{} - {}'.format(t1.strftime('%Y.%m.%d %H:%M'),t2.strftime('%Y.%m.%d %H:%M')))
for array in config.arrays:
    print('--- ' + array['Name'] + ' ---')
    #### download data ####
    SCNL  = DataFrame.from_dict(array['SCNL'])
    st    = utils.grab_data(SCNL['scnl'].tolist(),t1-config.latency,t2+config.latency+config.window_length,fill_value=0)
    st    = utils.add_coordinate_info(st,SCNL)
    array = utils.get_volcano_backazimuth(st,array)
    ########################

    #### check for enough data ####
    for tr in st:
        if np.sum(np.abs(tr.data))==0:
            st.remove(tr)
    if len(st)<config.min_chan:
        print('Too many blank traces. Skipping.')
        continue
    ########################

    #### check for gappy data ####
    for tr in st:
        if np.any([np.any(tr.data==0)]):
            st.remove(tr)
    if len(st)<config.min_chan:
        print('Too gappy. Skipping.')
        continue
    ########################

    #### preprocess data ####
    st.detrend('demean')
    st.taper(max_percentage=None,max_length=config.taper_val)
    st.filter('bandpass',freqmin=config.f1,freqmax=config.f2,corners=2,zerophase=True)
    st.trim(t1,t2+config.window_length+1)
    for tr in st:
        if tr.stats['sampling_rate']==100:
            tr.decimate(2)
        if tr.stats['sampling_rate']!=50:
            tr.resample(50.0)
        if tr.stats['sampling_rate']==50:
            tr.decimate(2)
    ########################


    velocity = []
    azimuth  = []
    mccm     = []
    t        = []
    for st_win in st.slide(window_length=config.window_length,step=config.overlap):
        try:
            for tr in st_win:
                tr.data = tr.data*array['digouti']
            vel, az, rms, cmax = utils.inversion(st_win)
            velocity.append(vel)
            azimuth.append(az)
            mccm.append(np.median(cmax))
            t_tmp=st_win[0].stats.starttime
            t.append(dates.date2num((st_win[0].stats.starttime+config.window_length/2.).datetime))
        except:
            continue
    t        = np.array(t)
    mccm     = np.array(mccm)
    velocity = np.array(velocity)
    azimuth  = np.array(azimuth)

    print('Making plot...')
    try:
        utils.web_folders(st,array,t2,config)
        utils.plot_results(t1,t2,t,st,mccm,velocity,azimuth,array,config)
    except:
        continue