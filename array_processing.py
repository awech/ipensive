#!/home/rtem/miniconda3/envs/ipensive/bin/python
# -*- coding: utf-8 -*-
"""
Created on 24-Apr-2018
Modified on 30-Mar-2022
@author: awech
"""

import os
import sys
sys.dont_write_bytecode = True  # don't write .pyc files (probably slightly faster without this, but more cluttered)
import numpy as np
import jinja2
from pandas import DataFrame
from obspy.core import UTCDateTime
from matplotlib import dates
import time
import config
import ipensive_utils as utils
import warnings
from matplotlib.cbook import MatplotlibDeprecationWarning
warnings.simplefilter('ignore',MatplotlibDeprecationWarning,append=True)
warnings.simplefilter('ignore',FutureWarning,append=True)


def process_array(array, network, T0):

	t1 = T0 - config.DURATION
	t2 = T0

	T1 = t1 - config.TAPER
	T2 = t2 + config.WINDOW_LENGTH + config.TAPER

	# get default params from config
	params={var:getattr(config,var) for var in dir(config) if var not in ['NETWORKS'] and '__' not in var}

	# update params with array-specific values
	for key in array.keys():
		if key in params.keys():
			params[key] = array[key]
	
	print('--- ' + array['Name'] + ' ---')
	time.sleep(params['EXTRA_PAUSE'])

	#### download data ####
	scnl  = DataFrame.from_dict(array['SCNL'])
	if len(scnl) < params['MIN_CHAN']:
		print('Not enough channels defined.')
		return
	st    = utils.grab_data(scnl['scnl'].tolist(), T1, T2,
							hostname=params['HOSTNAME'], port=params['PORT'], fill_value=0)
	st    = utils.add_coordinate_info(st,scnl)
	array = utils.get_volcano_backazimuth(st,array)
	########################

	#### check for enough data ####
	for tr in st:
		if np.sum(np.abs(tr.data))==0:
			st.remove(tr)
	if len(st) < params['MIN_CHAN']:
		print('Too many blank traces. Skipping.')
		return
	########################

	#### check for gappy data ####
	for tr in st:
		if np.any([np.any(tr.data==0)]):
			st.remove(tr)
	if len(st) < params['MIN_CHAN']:
		print('Too gappy. Skipping.')
		return
	########################

	#### preprocess data ####
	st.detrend('demean')
	for tr in st:
		if tr.stats['sampling_rate'] == 100:
			tr.decimate(2)
		if tr.stats['sampling_rate'] != 50:
			tr.resample(50.0)
		if tr.stats['sampling_rate'] == 50:
			tr.decimate(2)
	st.taper(max_percentage=None,max_length=params['TAPER'])
	st.filter('bandpass',freqmin=params['FREQMIN'],freqmax=params['FREQMAX'],corners=2,zerophase=True)
	st.trim(t1,t2+params['WINDOW_LENGTH'])
	########################

	velocity = []
	azimuth  = []
	mccm     = []
	t        = []
	rms      = []
	pressure = []
	for st_win in st.slide(window_length=params['WINDOW_LENGTH'],step=params['WINDOW_LENGTH']-params['OVERLAP']):
		try:
			for tr in st_win:
				tr.data = tr.data*array['digouti']
			vel, az, rms0, cmax, pk_press = utils.inversion(st_win)
			velocity.append(vel)
			azimuth.append(az)
			mccm.append(np.median(cmax))
			t_tmp=st_win[0].stats.starttime
			t.append(dates.date2num((st_win[0].stats.starttime+params['WINDOW_LENGTH']/2.).datetime))
			rms.append(rms0)
			pressure.append(pk_press)
		except:
			print('Something went wrong in the inversion...')
			continue
	t        = np.array(t)
	mccm     = np.array(mccm)
	velocity = np.array(velocity)
	azimuth  = np.array(azimuth)
	rms      = np.array(rms)
	pressure = np.array(pressure)

	try:
		print('Setting up web output folders')
		utils.web_folders(st,array,t2,network['Name'])
		print('Making plot...')
		utils.plot_results(t1, t2, t, st, mccm, velocity, azimuth, array, network['Name'])
	except:
		import traceback
		b=traceback.format_exc()
		message = ''.join('{}\n'.format(a) for a in b.splitlines())
		print('Something went wrong making the plot:')
		print(message)

	if 'OUT_VALVE_DIR' in dir(config):
		try:
			print('Writing csv file...')
			t=np.array([UTCDateTime(dates.num2date(ti)).strftime('%Y-%m-%d %H:%M:%S') for ti in t])
			name=st[0].stats.station
			utils.write_valve_file(t2, t, pressure, azimuth, velocity, mccm, rms, name)
		except:
			import traceback
			b=traceback.format_exc()
			message = ''.join('{}\n'.format(a) for a in b.splitlines())
			print('Something went wrong writing the csv file:')
			print(message)

	if 'OUT_ASCII_DIR' in dir(config):
		try:
			print('Writing csv file...')    
			t=np.array([UTCDateTime(dates.num2date(ti)).strftime('%Y-%m-%d %H:%M:%S') for ti in t])
			name=array['Name']
			utils.write_ascii_file(t2, t, pressure, azimuth, velocity, mccm, rms, name)
		except:
			import traceback
			b=traceback.format_exc()
			message = ''.join('{}\n'.format(a) for a in b.splitlines())
			print('Something went wrong writing the csv file:')
			print(message)

	return

if __name__ == '__main__':

	timer_start = time.time()
	if len(sys.argv) == 1:                              # No time given. Run from cron. Use current time
		
		T0=UTCDateTime.utcnow()                         # get current timestamp
		
		# Set up logging
		utils.write_to_log(T0.strftime('%Y-%m-%d'))
		
		# round down to the nearest 10-minute
		T0=UTCDateTime(T0.strftime('%Y-%m-%d %H:%M')[:-1]+'0')
		
		# Pause to allow for data latency to catch up
		time.sleep(config.LATENCY+config.WINDOW_LENGTH)
		timer_start = timer_start + config.LATENCY+config.WINDOW_LENGTH

	elif len(sys.argv) == 2:                            # time given as single string YYYYMMDDhhmm (eg. 201705130301)
		T0 = sys.argv[1] 
		try:
			T0 = UTCDateTime(T0)
			T0 = UTCDateTime(T0.strftime('%Y-%m-%d %H:%M')[:-1]+'0') # round down to the nearest 10-minute
		except:
			warnings.warn('Needs end-time argument. eg: array_processing.py 201701020205')
			sys.exit()

	print('Start time: {}'.format(UTCDateTime.utcnow()))
	print('Processing {} - {}'.format((T0 - config.DURATION).strftime('%Y.%m.%d %H:%M'),T0.strftime('%Y.%m.%d %H:%M')))

	all_nets=[]
	all_arrays={}
	for network in config.NETWORKS:
		all_nets.append(network['Name'])
		all_arrays[network['Name']] = []
		for array in network['ARRAYS']:
			all_arrays[network['Name']].append(array['Name'])
			timer_tmp = time.time()
			process_array(array, network, T0)
			print('{:.1f} seconds to process {}'.format(time.time()-timer_tmp, array['Name']))

	# Write out the new HTML file
	script_path = os.path.dirname(__file__)
	with open(os.path.join(script_path, 'index.template'), 'r') as f:
		template = jinja2.Template(f.read())
	html = template.render(networks = all_nets, arrays = all_arrays)
	with open(os.path.join(config.OUT_WEB_DIR, 'index.html'), 'w') as f:
		f.write(html)

	print('{:.1f} seconds to process all'.format(time.time() - timer_start))
	print('Finish time: {}'.format(UTCDateTime.utcnow()))
