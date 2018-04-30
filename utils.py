from os import path, mkdir
import sys
import numpy as np
from obspy.core import Stream
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

def add_coordinate_info(st,SCNL):
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
    rms = np.sqrt(np.mean(np.square(np.matmul(Dm3,sv)-dt.T)))

    return velocity, azimuth, rms, Cmax


def get_volcano_backazimuth(st,array):
    lon0=np.mean([tr.stats.coordinates.longitude for tr in st])
    lat0=np.mean([tr.stats.coordinates.latitude for tr in st])
    for volc in array['volcano']:
        tmp=gps2dist_azimuth(lat0,lon0,volc['v_lat'],volc['v_lon'])
        volc['back_azimuth']=tmp[1]
    return array

 
def grab_data(scnl,T1,T2,fill_value=0):
    # scnl = list of station names (eg. ['PS4A.EHZ.AV.--','PVV.EHZ.AV.--','PS1A.EHZ.AV.--'])
    # T1 and T2 are start/end obspy UTCDateTimes
    # fill_value can be 0 (default), 'latest', or 'interpolate'
    #
    # returns stream of traces with gaps accounted for
    #
    # print('{} - {}'.format(T1.strftime('%Y.%m.%d %H:%M:%S'),T2.strftime('%Y.%m.%d %H:%M:%S')))
    print('Grabbing data...')

    st=Stream()

    for sta in scnl:
        if sta.split('.')[2]=='MI':
            client = Client(config.winston_address_cnmi, config.winston_port_cnmi)
        else:
            client = Client(config.winston_address, config.winston_port)
        try:
            tr=client.get_waveforms(sta.split('.')[2], sta.split('.')[0],sta.split('.')[3],sta.split('.')[1], T1, T2, cleanup=True)
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
    st.trim(T1,T2,pad=0)
    print('Detrending data...')
    st.detrend('demean')
    return st


def web_folders(st,array,t2,config):
    from shutil import copyfile
    if not path.exists(config.out_dir):
        mkdir(config.out_dir)

    network=config.network
    if st[0].stats.network=='MI':
        network='CNMI'

    d0=config.out_dir+'/'+network
    if not path.exists(d0):
        mkdir(d0)
    d0=config.out_dir+'/'+network+'/'+array['Name']
    if not path.exists(d0):
        mkdir(d0)
    d0=config.out_dir+'/'+network+'/'+array['Name']+'/'+str(t2.year)
    if not path.exists(d0):
        mkdir(d0)
    d2=d0+'/'+'{:03d}'.format(t2.julday)
    if not path.exists(d2):
        mkdir(d2)
    copyfile(config.working_dir+'/index.html',config.out_dir+'/index.html')
    return


def plot_results(t1,t2,t,st,mccm,velocity,azimuth,array,config):

    ########## big plot ##########
    ##############################
    tvec     = np.linspace(dates.date2num(st[0].stats.starttime.datetime),dates.date2num(st[0].stats.endtime.datetime),len(st[0].data))
    T1=dates.date2num(t1.datetime)
    T2=dates.date2num(t2.datetime)
    cm='RdYlBu_r'
    cax=0.2,1   #colorbar/y-axis for mccm

    
    fig1=plt.figure(figsize=(8,10.5))
    axs1=plt.subplot(4,1,1)
    plt.title(array['Name']+' Infrasound Array')
    axs1.plot(tvec,st[0].data*array['digouti'],'k',linewidth=0.6)
    axs1.axis('tight')
    axs1.set_xlim(T1,T2)
    ymax=np.abs(list(axs1.get_ylim())).max()
    axs1.set_ylim(-ymax,ymax)
    # axs1.set_xticks([])
    axs1.xaxis_date()
    axs1.fmt_xdata = dates.DateFormatter('%HH:%MM')
    axs1.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
    axs1.set_xticklabels([])
    axs1.tick_params(direction='in',axis='x',top='on')
    axs1.set_ylabel('Pressure [Pa]')
    
    axs1=plt.subplot(4,1,2)
    sc=plt.scatter(t,mccm,c=mccm,edgecolors='k',lw=.3,cmap=cm)
    axs1.plot([T1,T2],[config.mcthresh,config.mcthresh],'--',color='gray')
    axs1.axis('tight')
    axs1.set_xlim(T1,T2)
    axs1.set_ylim(cax)
    sc.set_clim(cax)
    # axs1.set_xticks([])
    axs1.xaxis_date()
    axs1.fmt_xdata = dates.DateFormatter('%HH:%MM')
    axs1.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
    axs1.set_xticklabels([])
    axs1.tick_params(direction='in',axis='x',top='on')
    axs1.set_ylabel('MCCM')
    
    axs1=plt.subplot(4,1,3)
    rect=Rectangle((T1,0.25),T2-T1,0.45-0.25)
    pc = PatchCollection([rect], facecolor='gray', alpha=0.25,edgecolor=None)
    plt.gca().add_collection(pc)
    sc=axs1.scatter(t,velocity,c=mccm,edgecolors='k',lw=.3,cmap=cm)
    axs1.set_ylim(.15,.6)
    axs1.set_xlim(T1,T2)
    sc.set_clim(cax)
    # axs1.set_xticks([])
    axs1.xaxis_date()
    axs1.fmt_xdata = dates.DateFormatter('%HH:%MM')
    axs1.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
    axs1.set_xticklabels([])
    axs1.tick_params(direction='in',axis='x',top='on')
    axs1.set_ylabel('Trace Velocity\n [km/s]')
    
    axs1=plt.subplot(4,1,4)
    for volc in array['volcano']:
        axs1.plot([T1,T2],[volc['back_azimuth'],volc['back_azimuth']],'--',color='gray',zorder=-1)
        axs1.text(t[1],volc['back_azimuth']-6,volc['name'],bbox={'facecolor':'white','edgecolor':'white','pad':0},fontsize=8,style='italic',zorder=10)
    sc=axs1.scatter(t,azimuth,c=mccm,edgecolors='k',lw=.3,cmap=cm,zorder=1000)

    axs1.set_ylim(0,360)
    axs1.set_xlim(T1,T2)
    sc.set_clim(cax)
    axs1.set_ylabel('Back-azimuth\n [deg]')
    
    axs1.xaxis_date()
    axs1.tick_params(axis='x',labelbottom='on')
    axs1.fmt_xdata = dates.DateFormatter('%HH:%MM')
    axs1.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
    axs1.set_xlabel('UTC Time ['+t1.strftime('%Y-%b-%d')+']')
    
    plt.subplots_adjust(left=0.1,right=.9,top=0.97,bottom=0.05,hspace=0.1)
    axs1=plt.subplot(4,1,2)
    ctop=axs1.get_position().y1
    axs1=plt.subplot(4,1,4)
    cbot=axs1.get_position().y0
    cbaxes=fig1.add_axes([.91,cbot,.02,ctop-cbot])
    hc=plt.colorbar(sc,cax=cbaxes)
    hc.set_label('MCCM')
    
    network=config.network
    if st[0].stats.network=='MI':
        network='CNMI'
    d0=config.out_dir+'/'+network+'/'+array['Name']+'/'+str(t2.year)
    d2=d0+'/'+'{:03d}'.format(t2.julday)
    filename=d2+'/'+array['Name']+'_'+t2.strftime('%Y%m%d-%H%M')+'.png'
    plt.savefig(filename,dpi=72,format='png')
    plt.close('all')



    ######### small plot #########
    ##############################
    fig1=plt.figure(figsize=(2.1,2.75))
    axs1=plt.subplot(4,1,1)
    axs1.plot(tvec,st[0].data*array['digouti'],'k',linewidth=0.1)
    axs1.axis('tight')
    ymax=np.abs(list(axs1.get_ylim())).max()
    axs1.set_ylim(-ymax,ymax)
    axs1.set_xlim(T1,T2)
    axs1.set_xticks([])
    axs1.set_yticks([])
    
    axs1=plt.subplot(4,1,2)
    sc=plt.scatter(t,mccm,s=8*np.ones_like(t),c=mccm,edgecolors='k',lw=.1,cmap=cm)
    axs1.plot([T1,T2],[config.mcthresh,config.mcthresh],'--',color='gray',linewidth=1)
    axs1.axis('tight')
    axs1.set_xlim(T1,T2)
    axs1.set_ylim(cax)
    sc.set_clim(cax)
    axs1.set_xticks([])
    axs1.set_yticks([])
    
    axs1=plt.subplot(4,1,3)
    rect=Rectangle((T1,0.25),T2-T1,0.45-0.25)
    pc = PatchCollection([rect], facecolor='gray', alpha=0.25,edgecolor=None)
    plt.gca().add_collection(pc)
    sc=axs1.scatter(t,velocity,s=8*np.ones_like(t),c=mccm,edgecolors='k',lw=.1,cmap=cm)
    axs1.set_ylim(.15,.6)
    axs1.set_xlim(T1,T2)
    sc.set_clim(cax)
    axs1.set_xticks([])
    axs1.set_yticks([])

    
    axs1=plt.subplot(4,1,4)
    sc=axs1.scatter(t,azimuth,s=8*np.ones_like(t),c=mccm,edgecolors='k',lw=.3,cmap=cm)
    axs1.set_ylim(0,360)
    axs1.set_xlim(T1,T2)
    sc.set_clim(cax)
    axs1.set_xticks([])
    axs1.set_yticks([])

    plt.subplots_adjust(left=0,right=1,top=0.99,bottom=0.01,hspace=0.03)
    filename=d2+'/'+array['Name']+'_'+t2.strftime('%Y%m%d-%H%M')+'_thumb.png'
    plt.savefig(filename,format='png',pad_inches=0,dpi=72)
    plt.close('all')
