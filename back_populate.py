#!/Users/awech-local/opt/miniconda3/envs/ipensive/bin/python
# -*- coding: utf-8 -*-
"""
Created on 30-Mar-2022
@author: awech
"""

import os
import pandas as pd
from obspy import UTCDateTime
import config
from array_processing import process_array

T1='2022-03-30 12:00'
T2='2022-03-30 15:00'
OVERWRITE=False

# ARRAYS=['Dillingham','Kenai','Sand Point','Akutan','Cleveland','Adak','Whittier']
# specify ARRAYS if you want to just process specific arrays

def run_backpopulate():

	t1 = UTCDateTime(T1)+config.DURATION
	for t in pd.date_range(T2, t1.strftime('%Y%m%d%H%M'), freq='-10min'):
		for network in config.NETWORKS:
			for array in network['ARRAYS']:

				# set up filename
				year='{:.0f}'.format(t.year)
				day='{:03.0f}'.format(t.dayofyear)
				time=t.strftime('%Y%m%d-%H%M')
				file='{}/{}/{}/{}/{}/{}_{}.png'.format(config.OUT_WEB_DIR,
													   network['Name'],
													   array['Name'],
													   year,
													   day,
													   array['Name'],
													   time)

				# check if you should process this time window
				if not os.path.exists(file) or (os.path.exists(file) and OVERWRITE):
					
					# check if you should process this array
					if 'ARRAYS' in dir():
						if array['Name'] in ARRAYS:
							process_array(array, network, UTCDateTime(t))
					else:
						process_array(array, network, UTCDateTime(t))


if __name__ == '__main__':

	run_backpopulate()
