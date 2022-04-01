import os
import sys
import numpy as np
import pandas as pd
from obspy.core import Stream, UTCDateTime
from obspy.core.util import AttribDict
from obspy.clients.earthworm import Client
from obspy.geodetics.base import gps2dist_azimuth
from scipy.signal import correlate
import config


####### plotting imports #######
import matplotlib as m
m.use('Agg')
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib import dates
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
fonts=10
rcParams.update({'font.size': fonts})
################################


def write_to_log(day):

	if 'LOGS_DIR' in dir(config) and len(config.LOGS_DIR)>0:
		LOGS_DIR = config.LOGS_DIR
	else:
		LOGS_DIR = os.path.dirname(__file__) + '/logs'

	year=UTCDateTime(day).strftime('%Y')
	month=UTCDateTime(day).strftime('%Y-%m')
	if not os.path.exists(LOGS_DIR+'/'+year):
		os.mkdir(LOGS_DIR+'/'+year)
	if not os.path.exists(LOGS_DIR+'/'+year+'/'+month):
		os.mkdir(LOGS_DIR+'/'+year+'/'+month)
	file=LOGS_DIR+'/'+year+'/'+month+'/'+UTCDateTime(day).strftime('%Y-%m-%d')+'.log'
	os.system('touch {}'.format(file))
	f=open(file,'a')
	sys.stdout=sys.stderr=f
	return


def add_coordinate_info(st, SCNL):
	#### compare remaining stations with lat/lon station info in config file
	#### to attach lat/lon info with each corresponding trace
	for tr in st:
		if tr.stats.location=='':
			tr.stats.location='--'
		tmp_scnl='{}.{}.{}.{}'.format(tr.stats.station,
									  tr.stats.channel,
									  tr.stats.network,
									  tr.stats.location)
		tmp_lat=SCNL[SCNL['scnl']==tmp_scnl].sta_lat.values[0]
		tmp_lon=SCNL[SCNL['scnl']==tmp_scnl].sta_lon.values[0]
		tr.stats.coordinates=AttribDict({
								'latitude': tmp_lat,
								'longitude': tmp_lon,
								'elevation': 0.0})
	return st


def setup_coordinate_system(st):
	R = 6372.7976   # radius of the earth
	lons  = np.array([tr.stats.coordinates.longitude for tr in st])
	lats  = np.array([tr.stats.coordinates.latitude for tr in st])
	lon0  = lons.mean()*np.pi/180.0
	lat0  = lats.mean()*np.pi/180.0
	yx    = R*np.array([ lats*np.pi/180.0-lat0, (lons*np.pi/180.0-lon0)*np.cos(lat0) ]).T
	intsd = np.zeros([len(lons),len(lons)])
	ints_az= np.zeros([len(lons),len(lons)])
	for ii in range(len(st[:-1])):
		for jj in range(ii+1,len(st)):
			# intsd[i,j]=np.sqrt(np.square(yx[j][0]-yx[i][0])+np.square(yx[j][1]-yx[i][1]))
			tmp=gps2dist_azimuth(lats[ii],lons[ii],lats[jj],lons[jj])
			intsd[ii,jj]=tmp[0]
			ints_az[ii,jj]=tmp[1]

	return yx, intsd, ints_az


def align_stack_stream(st, LAGS):

	st_temp=st.copy()
	group_streams = Stream()
	T1 = st_temp[0].copy().stats.starttime
	T2 = st_temp[0].copy().stats.endtime
	for i, tr in enumerate(st_temp):
		tr = tr.copy().trim(
			tr.stats.starttime - LAGS[i],
			tr.stats.endtime - LAGS[i],
			pad=True,
			fill_value=0,
		)
		tr.trim(tr.stats.starttime + 1, tr.stats.endtime - 1, pad=True, fill_value=0)
		tr.stats.starttime = T1
		group_streams += tr

	ST = group_streams[0].copy()
	for tr in group_streams[1:]:
		ST.data = ST.data + tr.data
	ST.data = ST.data / len(st_temp)
	ST.trim(T1, T2)
	return ST


def inversion(st):
	## inversion originally written by M. Haney in matlab
	## modified and converted to Python by A. Wech
	lags   = np.array([])
	Cmax   = np.array([])

	mlag = st[0].stats.npts
	tC   = np.linspace(-mlag,mlag,2*mlag-1)/st[0].stats.sampling_rate
	for ii in range(len(st[:-1])):
		for jj in range(ii+1,len(st)):
			scale=np.linalg.norm(st[ii].data)*np.linalg.norm(st[jj].data)
			cc=correlate(st[ii],st[jj],mode='full')/float(scale)
			Cmax = np.append(Cmax,cc.max())
			lags = np.append(lags,tC[cc.argmax()])

	# get interstation distance and azimuth vectors
	yx, intsd, ints_az = setup_coordinate_system(st)
	ds = intsd[np.triu_indices(len(st),1)]
	az = ints_az[np.triu_indices(len(st),1)]

	dt    = lags
	Dm3   = np.array([ds*np.cos(az*(np.pi/180.0)) , ds*np.sin(az*(np.pi/180.0))]).T
	Dm3   = Dm3/1000.0  # convert to kilometers

	# generalized inverse of slowness matrix
	Gmi = np.linalg.inv(np.matmul(Dm3.T,Dm3))
	# slowness - least squares
	sv = np.matmul(np.matmul(Gmi,Dm3.T),dt.T)
	# velocity from slowness
	velocity = 1/np.sqrt(np.square(sv[0])+np.square(sv[1]))
	# cosine and sine for backazimuth
	caz3 = velocity*sv[0]
	saz3 = velocity*sv[1]
	# 180 degree resolved backazimuth to source
	azimuth = np.arctan2(saz3,caz3)*(180/np.pi)
	if azimuth<0:
		azimuth=azimuth+360
	# rms
	# rms = np.sqrt(np.mean(np.square(np.matmul(Dm3,sv)-dt.T)))

	Dm3_new=np.array([ds*np.cos(az*(np.pi/180.0)) , ds*np.sin(az*(np.pi/180.0))]).T/1000
	sv_new=np.array([np.cos(azimuth*np.pi/180)/velocity, np.sin(azimuth*np.pi/180)/velocity])
	lags_new=np.matmul(Dm3_new,sv_new)
	rms = np.sqrt(np.mean(np.square(lags_new-dt.T)))

	LAGS=np.hstack((0,lags_new[:len(st)-1]))
	ST_stack=align_stack_stream(st,LAGS)
	pk_pressure=np.max(np.abs(ST_stack.data))

	return velocity, azimuth, rms, Cmax, pk_pressure


def get_volcano_backazimuth(st, array):
	lon0=np.mean([tr.stats.coordinates.longitude for tr in st])
	lat0=np.mean([tr.stats.coordinates.latitude for tr in st])
	for volc in array['volcano']:
		if 'back_azimuth' not in volc:
			tmp=gps2dist_azimuth(lat0,lon0,volc['v_lat'],volc['v_lon'])
			volc['back_azimuth']=tmp[1]
	return array

 
def grab_data(scnl, T1, T2, hostname, port, fill_value=0):
	# scnl = list of station names (eg. ['PS4A.EHZ.AV.--','PVV.EHZ.AV.--','PS1A.EHZ.AV.--'])
	# T1 and T2 are start/end obspy UTCDateTimes
	# fill_value can be 0 (default), 'latest', or 'interpolate'
	#
	# returns stream of traces with gaps accounted for
	#
	# print('{} - {}'.format(T1.strftime('%Y.%m.%d %H:%M:%S'),T2.strftime('%Y.%m.%d %H:%M:%S')))
	print('Grabbing data...')

	st=Stream()
	client = Client(hostname, int(port))

	for sta in scnl:
		
		try:
			tr=client.get_waveforms(sta.split('.')[2], sta.split('.')[0],sta.split('.')[3],sta.split('.')[1],
									T1, T2, cleanup=True)
			if len(tr)>1:
				if fill_value==0 or fill_value==None:
					tr.detrend('demean')
					tr.taper(max_percentage=0.01)
				for sub_trace in tr:
					# deal with error when sub-traces have different dtypes
					if sub_trace.data.dtype.name != 'int32':
						sub_trace.data=sub_trace.data.astype('int32')
					if sub_trace.data.dtype!=np.dtype('int32'):
						sub_trace.data=sub_trace.data.astype('int32')
					# deal with rare error when sub-traces have different sample rates
					if sub_trace.stats.sampling_rate!=np.round(sub_trace.stats.sampling_rate):
						sub_trace.stats.sampling_rate=np.round(sub_trace.stats.sampling_rate)
				print('Merging gappy data...')
				tr.merge(fill_value=fill_value)

			# deal where trace length is smaller than expected window length
			if tr[0].stats.endtime - tr[0].stats.starttime < T2 - T1:
				tr.detrend('demean')
				tr.taper(max_percentage=0.01)
		except:
			tr=Stream()
		# if no data, create a blank trace for that channel
		if not tr:
			from obspy import Trace
			from numpy import zeros
			tr=Trace()
			tr.stats['station']=sta.split('.')[0]
			tr.stats['channel']=sta.split('.')[1]
			tr.stats['network']=sta.split('.')[2]
			tr.stats['location']=sta.split('.')[3]
			tr.stats['sampling_rate']=100
			tr.stats['starttime']=T1
			tr.data=zeros(int((T2-T1)*tr.stats['sampling_rate']),dtype='int32')
		st+=tr
	st.trim(T1,T2,pad=True, fill_value=0)
	print('Detrending data...')
	st.detrend('demean')
	return st


def web_folders(st, array, t2, network):
	from shutil import copyfile
	if not os.path.exists(config.OUT_WEB_DIR):
		os.mkdir(config.OUT_WEB_DIR)

	d0=config.OUT_WEB_DIR+'/'+network
	if not os.path.exists(d0):
		os.mkdir(d0)
	d0=config.OUT_WEB_DIR+'/'+network+'/'+array['Name']
	if not os.path.exists(d0):
		os.mkdir(d0)
	d0=config.OUT_WEB_DIR+'/'+network+'/'+array['Name']+'/'+str(t2.year)
	if not os.path.exists(d0):
		os.mkdir(d0)
	d2=d0+'/'+'{:03d}'.format(t2.julday)
	if not os.path.exists(d2):
		os.mkdir(d2)

	# copyfile('index.html',config.OUT_WEB_DIR+'/index.html')
	return


def write_ascii_file(t2, t, pressure, azimuth, velocity, mccm, rms, name):

	name=name.replace(' ','_')

	t1=t2-config.DURATION

	d0=config.OUT_ASCII_DIR+'/'+name
	if not os.path.exists(d0):
		os.mkdir(d0)

	subfolder=d0+'/{}'.format(t1.strftime('%Y-%m'))
	if not os.path.exists(subfolder):
		os.mkdir(subfolder)

	
	filename=subfolder+'/'+name+'_'+t1.strftime('%Y-%m-%d')+'.txt'

	azimuth[azimuth<0]+=360

	tmp=pd.DataFrame({'Time':t,
				 'Array':name,
				 'Azimuth':azimuth,
				 'Velocity':velocity,
				 'MCCM':mccm,
				 'Pressure':pressure,
				 'rms':rms})
	tmp['Time']=pd.to_datetime(tmp['Time'])
	tmp = tmp[tmp['Time']<=t2.strftime('%Y-%m-%d %H:%M:%S')]
	tmp['Velocity']=1000*tmp['Velocity']

	if os.path.exists(filename):
		df = pd.read_csv(filename, sep='\t', parse_dates=['Time'])
		df = df[(df['Time'] <= t1.strftime('%Y-%m-%d %H:%M:%S')) | (df['Time'] > t2.strftime('%Y-%m-%d %H:%M:%S'))]
		df = pd.concat([df,tmp])
		df = df.sort_values('Time')
	else:
		df = tmp

	df = df.round({'Azimuth':1,'Velocity':1,'MCCM':2,'Pressure':3,'rms':1})

	df.to_csv(filename,index=False,header=True,sep='\t')
	return


def write_valve_file(t2, t, pressure, azimuth, velocity, mccm, rms, name):
	from pandas import DataFrame

	A=DataFrame({'TIMESTAMP':t,
				 'CHANNEL':name,
				 'Azimuth':azimuth,
				 'Velocity':velocity,
				 'MCCM':mccm,
				 'Pressure':pressure,
				 'rms':rms})

	A=A[['TIMESTAMP','CHANNEL','Azimuth','Velocity','MCCM','Pressure','rms']]

	A['Velocity']=1000*A['Velocity']

	filename=config.out_valve_dir+'/'+name+'_'+t2.strftime('%Y%m%d-%H%M')+'.txt'

	A.to_csv(filename,index=False,header=True,sep=',',float_format='%.3f')


def plot_results(t1, t2, t, st, mccm, velocity, azimuth, array, network):


	# get default params from config
	params_tmp={var:getattr(config,var) for var in dir(config) if var not in ['NETWORKS'] and '__' not in var}

	# update params with array-specific values
	for key in array.keys():
		if key in params_tmp.keys():
			params_tmp[key] = array[key]
	

	########## big plot ##########
	##############################
	tvec     = np.linspace(dates.date2num(st[0].stats.starttime.datetime),dates.date2num(st[0].stats.endtime.datetime),len(st[0].data))
	T1=dates.date2num(t1.datetime)
	T2=dates.date2num(t2.datetime)
	cm='RdYlBu_r'
	cax=0.2,1   #colorbar/y-axis for mccm

	
	fig1=plt.figure(figsize=(8,10.5))
	axs1=plt.subplot(4,1,1)
	plt.title(array['Name']+' '+params_tmp['ARRAY_LABEL']+ ' Array')
	axs1.plot(tvec,st[0].data*array['digouti'],'k',linewidth=0.6)
	axs1.axis('tight')
	axs1.set_xlim(T1,T2)
	ymax=np.abs(list(axs1.get_ylim())).max()
	axs1.set_ylim(-ymax,ymax)
	axs1.xaxis_date()
	axs1.fmt_xdata = dates.DateFormatter('%HH:%MM')
	axs1.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
	axs1.set_xticklabels([])
	axs1.tick_params(direction='in',axis='x',top='on')
	axs1.set_ylabel('Pressure [Pa]')
	axs1b = axs1.twinx()
	axs1b.set_yticks([])
	axs1b.set_ylabel('{:.1f} - {:.1f} Hz'.format(params_tmp['FREQMIN'], params_tmp['FREQMAX']), labelpad=6)
	
	axs2=plt.subplot(4,1,2)
	sc=plt.scatter(t,mccm,c=mccm,edgecolors='k',lw=.3,cmap=cm)
	axs2.plot([T1,T2],[params_tmp['MCTHRESH'], params_tmp['MCTHRESH']],'--',color='gray')
	axs2.axis('tight')
	axs2.set_xlim(T1,T2)
	axs2.set_ylim(cax)
	sc.set_clim(cax)
	axs2.xaxis_date()
	axs2.fmt_xdata = dates.DateFormatter('%HH:%MM')
	axs2.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
	axs2.set_xticklabels([])
	axs2.tick_params(direction='in',axis='x',top='on')
	axs2.set_ylabel('MCCM')
	
	axs3=plt.subplot(4,1,3)
	rect=Rectangle((T1,params_tmp['VEL_MIN']),T2-T1,params_tmp['VEL_MAX'] - params_tmp['VEL_MIN'])
	pc = PatchCollection([rect], facecolor='gray', alpha=0.25,edgecolor=None)
	plt.gca().add_collection(pc)
	sc=axs3.scatter(t,velocity,c=mccm,edgecolors='k',lw=.3,cmap=cm)
	if params_tmp['ARRAY_LABEL'] == 'Hydroacoustic':
		axs3.set_ylim(1.2,1.8)
	else:
		axs3.set_ylim(.15,.6)
	axs3.set_xlim(T1,T2)
	sc.set_clim(cax)
	axs3.xaxis_date()
	axs3.fmt_xdata = dates.DateFormatter('%HH:%MM')
	axs3.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
	axs3.set_xticklabels([])
	axs3.tick_params(direction='in',axis='x',top='on')
	axs3.set_ylabel('Trace Velocity\n [km/s]')
	
	axs4=plt.subplot(4,1,4)

	if params_tmp['AZ_MAX'] < params_tmp['AZ_MIN']:
		params_tmp['AZ_MIN'] = params_tmp['AZ_MIN']-360
		for volc in array['volcano']:
			if volc['back_azimuth']>180:
				axs4.plot([T1,T2],[volc['back_azimuth']-360,volc['back_azimuth']-360],'--',color='gray',zorder=-1)
				axs4.text(t[1],volc['back_azimuth']-360,volc['name'],bbox={'facecolor':'white','edgecolor':'white','pad':0},fontsize=8,verticalalignment='center',style='italic',zorder=10)
			else:
				axs4.plot([T1,T2],[volc['back_azimuth'],volc['back_azimuth']],'--',color='gray',zorder=-1)
				axs4.text(t[1],volc['back_azimuth'],volc['name'],bbox={'facecolor':'white','edgecolor':'white','pad':0},fontsize=8,verticalalignment='center',style='italic',zorder=10)
		azimuth[azimuth>180]+=-360
	else:
		for volc in array['volcano']:
			axs4.plot([T1,T2],[volc['back_azimuth'],volc['back_azimuth']],'--',color='gray',zorder=-1)
			axs4.text(t[1],volc['back_azimuth'],volc['name'],bbox={'facecolor':'white','edgecolor':'white','pad':0},fontsize=8,verticalalignment='center',style='italic',zorder=10)
	sc=axs4.scatter(t,azimuth,c=mccm,edgecolors='k',lw=.3,cmap=cm,zorder=1000)
	axs4.set_ylim(params_tmp['AZ_MIN'],params_tmp['AZ_MAX'])
	axs4.set_xlim(T1,T2)
	sc.set_clim(cax)
	axs4.set_ylabel('Back-azimuth\n [deg]')
	
	axs4.xaxis_date()
	axs4.tick_params(axis='x',labelbottom='on')
	axs4.fmt_xdata = dates.DateFormatter('%HH:%MM')
	axs4.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
	axs4.set_xlabel('UTC Time ['+t1.strftime('%Y-%b-%d')+']')
	
	plt.subplots_adjust(left=0.1,right=.9,top=0.97,bottom=0.05,hspace=0.1)
	ctop=axs2.get_position().y1
	cbot=axs4.get_position().y0
	cbaxes=fig1.add_axes([.91,cbot,.02,ctop-cbot])
	hc=plt.colorbar(sc,cax=cbaxes)
	hc.set_label('MCCM')
	
	print(array['Name'])
	d0=config.OUT_WEB_DIR+'/'+network+'/'+array['Name']+'/'+str(t2.year)
	d2=d0+'/'+'{:03d}'.format(t2.julday)
	filename=d2+'/'+array['Name']+'_'+t2.strftime('%Y%m%d-%H%M')+'.png'
	plt.savefig(filename,dpi=72,format='png')
	plt.close('all')



	######### small plot #########
	##############################
	fig1=plt.figure(figsize=(2.1,2.75))
	ax_1=plt.subplot(4,1,1)
	ax_1.plot(tvec,st[0].data*array['digouti'],'k',linewidth=0.1)
	ax_1.axis('tight')
	ymax=np.abs(list(ax_1.get_ylim())).max()
	ax_1.set_ylim(-ymax,ymax)
	ax_1.set_xlim(T1,T2)
	ax_1.set_xticks([])
	ax_1.set_yticks([])
	
	ax_2=plt.subplot(4,1,2)
	sc=plt.scatter(t,mccm,s=8*np.ones_like(t),c=mccm,edgecolors='k',lw=.1,cmap=cm)
	ax_2.plot([T1,T2],[params_tmp['MCTHRESH'], params_tmp['MCTHRESH']],'--',color='gray',linewidth=1)
	ax_2.axis('tight')
	ax_2.set_xlim(T1,T2)
	ax_2.set_ylim(cax)
	sc.set_clim(cax)
	ax_2.set_xticks([])
	ax_2.set_yticks([])
	
	ax_3=plt.subplot(4,1,3)
	rect=Rectangle((T1,0.25),T2-T1,params_tmp['VEL_MAX'] - params_tmp['VEL_MIN'])
	pc = PatchCollection([rect], facecolor='gray', alpha=0.25,edgecolor=None)
	plt.gca().add_collection(pc)
	sc=ax_3.scatter(t,velocity,s=8*np.ones_like(t),c=mccm,edgecolors='k',lw=.1,cmap=cm)
	ax_3.set_ylim(.15,.6)
	ax_3.set_xlim(T1,T2)
	sc.set_clim(cax)
	ax_3.set_xticks([])
	ax_3.set_yticks([])
	
	ax_4=plt.subplot(4,1,4)
	if params_tmp['AZ_MAX'] < params_tmp['AZ_MIN']:
		params_tmp['AZ_MIN'] = params_tmp['AZ_MIN']-360
		azimuth[azimuth>180]+=-360

	sc=ax_4.scatter(t,azimuth,s=8*np.ones_like(t),c=mccm,edgecolors='k',lw=.3,cmap=cm)
	ax_4.set_ylim(params_tmp['AZ_MIN'],params_tmp['AZ_MAX'])
	ax_4.set_xlim(T1,T2)
	sc.set_clim(cax)
	ax_4.set_xticks([])
	ax_4.set_yticks([])

	plt.subplots_adjust(left=0,right=1,top=0.99,bottom=0.01,hspace=0.03)
	filename=d2+'/'+array['Name']+'_'+t2.strftime('%Y%m%d-%H%M')+'_thumb.png'
	plt.savefig(filename,format='png',pad_inches=0,dpi=72)
	plt.close('all')
