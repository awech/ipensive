import os
import sys
import numpy as np
import pandas as pd
from obspy import Stream, UTCDateTime, read_inventory
from obspy.clients.earthworm import Client
from obspy.clients.fdsn import Client as FDSNClient
from obspy.geodetics.base import gps2dist_azimuth
from scipy.signal import correlate
import yaml

####### plotting imports #######
import matplotlib as m
m.use('Agg')
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib import dates
fonts=10
rcParams.update({'font.size': fonts})
################################


def load_config(config_file):

    with open(config_file, "r") as file:
        config = yaml.safe_load(file)

    all_nets = list(config["NETWORKS"].keys())
    array_list = []
    for net in all_nets:
        for array in config["NETWORKS"][net]:
            config[array]["NETWORK_NAME"] = net
            config[array]["ARRAY_NAME"] = array
            array_list.append(array)

    config["network_list"] = all_nets
    config["array_list"] = array_list

    return config


def write_to_log(day, config):

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


def check_inventory(tr, inv):
    inv_test = inv.select(
        network=tr.stats.network,
        station=tr.stats.station,
        location=tr.stats.location,
        channel=tr.stats.channel,
        starttime=tr.stats.starttime,
        endtime=tr.stats.starttime,
    )
    value = True if len(inv_test) > 0 else False
    return value


def check_FDSN(tr, client):
    value = True
    try:
        inventory = client.get_stations(
            network=tr.stats.network,
            station=tr.stats.station,
            location=tr.stats.location,
            channel=tr.stats.channel,
            starttime=tr.stats.starttime,
            endtime=tr.stats.starttime,
            level="response",
        )
    except Exception as err:
        if "No data available for request." in err.args[0]:
            value = False
    return value


def add_metadata(st, xml_file_name):
    import warnings

    warnings.simplefilter("ignore", UserWarning, append=True)
    inventory = read_inventory(xml_file_name)

    for tr in st:
        print(f"Getting metadata for {tr.id}")
        if check_inventory(tr, inventory):
            inv = inventory.select(
                network=tr.stats.network,
                station=tr.stats.station,
                location=tr.stats.location,
                channel=tr.stats.channel,
                starttime=tr.stats.starttime,
                endtime=tr.stats.endtime,
            )
        else:
            print(
                f"No station response info in stationXML file. Getting station response for {tr.id} from IRIS"
            )
            client = FDSNClient("IRIS")
            if check_FDSN(tr, client):
                inv = client.get_stations(
                    network=tr.stats.network,
                    station=tr.stats.station,
                    location=tr.stats.location,
                    channel=tr.stats.channel,
                    starttime=tr.stats.starttime,
                    endtime=tr.stats.endtime,
                    level="response",
                )
            else:
                print(f"No data available for request. Removing {tr.id}")
                st.remove(tr)
                continue
        tr.stats.coordinates = inv.get_coordinates(tr.id, tr.stats.starttime)
        tr.inventory = inv

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


def get_volcano_backazimuth(st, config, params):
    # def get_volcano_backazimuth(st, array):
    lon0=np.mean([tr.stats.coordinates.longitude for tr in st])
    lat0=np.mean([tr.stats.coordinates.latitude for tr in st])
    DF = pd.read_csv(config["TARGETS_FILE"])
    for target in params["TARGETS"]:
        if type(target) is str:
            df = DF[DF["Target"] == target]
            _, baz, _ = gps2dist_azimuth(lat0, lon0, df.iloc[0]["Latitude"], df.iloc[0]["Longitude"])
            params[target] = baz
    return params


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

    if hostname == 'IRIS':
        client = FDSNClient('IRIS')
    else:
        client = Client(hostname, int(port))

    for sta in scnl:
        
        try:
            if hostname == 'IRIS':
                tr=client.get_waveforms(sta.split('.')[2], sta.split('.')[0],sta.split('.')[3],sta.split('.')[1],
                                    T1, T2)
            else:
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


def web_folders(t2, config, params):

    if not os.path.exists(config["OUT_WEB_DIR"]):
        os.mkdir(config["OUT_WEB_DIR"])

    d0 = config["OUT_WEB_DIR"] + "/" + params["NETWORK_NAME"]
    if not os.path.exists(d0):
        os.mkdir(d0)
    d0 = config["OUT_WEB_DIR"] + "/" + params["NETWORK_NAME"] + "/" + params["ARRAY_NAME"]
    if not os.path.exists(d0):
        os.mkdir(d0)
    d0 = (
        config["OUT_WEB_DIR"] + "/" + params["NETWORK_NAME"] + "/" + params["ARRAY_NAME"] + "/" + str(t2.year)
    )
    if not os.path.exists(d0):
        os.mkdir(d0)
    d2 = d0 + "/" + "{:03d}".format(t2.julday)
    if not os.path.exists(d2):
        os.mkdir(d2)

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


def plot_results(t1, t2, t, st, mccm, velocity, azimuth, config, params):

    d0=config["OUT_WEB_DIR"]+'/'+params["NETWORK_NAME"]+'/'+params["ARRAY_NAME"]+'/'+str(t2.year)
    d2=d0+'/'+'{:03d}'.format(t2.julday)
    
    tvec = np.linspace(
        dates.date2num(st[0].stats.starttime.datetime),
        dates.date2num(st[0].stats.endtime.datetime),
        len(st[0].data),
    )
    T1 = dates.date2num(t1.datetime)
    T2 = dates.date2num(t2.datetime)
    cm = "RdYlBu_r"
    cax = 0.2, 1  # colorbar/y-axis for mccm

    size = (8, 10.5)
    trace_lw = 0.6
    scatter_lw = 0.3
    hline_lw = 1
    for plot_size in ["big", "small"]:
        if plot_size == "small":
            size = (2.1, 2.75)
            trace_lw = 0.1
            scatter_lw = 0.1
            hline_lw = 0.5

        fig, ax = plt.subplot_mosaic(
            [["wave"], ["cc"], ["vel"], ["baz"]],
            figsize=size
        )

        ############ Plot Waveforms ##############
        ##########################################
        ax["wave"].set_title(params["ARRAY_NAME"] + " " + params["ARRAY_LABEL"] + " Array")
        ax["wave"].plot(tvec, st[0].data, "k", linewidth=trace_lw)
        ax["wave"].axis("tight")
        ax["wave"].set_xlim(T1, T2)
        ymax = np.abs(list(ax["wave"].get_ylim())).max()
        ax["wave"].set_ylim(-ymax, ymax)
        if plot_size == "big":
            ax["wave"].xaxis_date()
            ax["wave"].fmt_xdata = dates.DateFormatter("%HH:%MM")
            ax["wave"].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
            ax["wave"].set_xticklabels([])
            ax["wave"].tick_params(direction="in", axis="x", top="on")
            ax["wave"].set_ylabel("Pressure [Pa]")
            ax["wave2"] = ax["wave"].twinx()
            ax["wave2"].set_yticks([])
            ax["wave2"].set_ylabel(
                f"{params["FREQMIN"]:.1f} - {params["FREQMAX"]:.1f} Hz",
                labelpad=6,
            )
        else:
            ax["wave"].set_xticks([])
            ax["wave"].set_yticks([])
        ##########################################


        ############ Plot cc values ##############
        ##########################################
        if plot_size == "big":
            sc = ax["cc"].scatter(t, mccm, c=mccm, edgecolors="k", lw=scatter_lw, cmap=cm)
        else:
            sc = ax["cc"].scatter(t, mccm, s=8*np.ones_like(t), c=mccm, edgecolors="k", lw=scatter_lw, cmap=cm)
        ax["cc"].axhline(params["MCTHRESH"], ls="--", lw=hline_lw, color="gray")
        ax["cc"].axis("tight")
        ax["cc"].set_xlim(T1, T2)
        ax["cc"].set_ylim(cax)
        sc.set_clim(cax)
        if plot_size == "big":
            ax["cc"].xaxis_date()
            ax["cc"].fmt_xdata = dates.DateFormatter("%HH:%MM")
            ax["cc"].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
            ax["cc"].set_xticklabels([])
            ax["cc"].tick_params(direction="in", axis="x", top="on")
            ax["cc"].set_ylabel(r"$M_{d}CCM$")
        else:
            ax["cc"].set_xticks([])
            ax["cc"].set_yticks([])
        ##########################################

        ########## Plot Trace Velocities #########
        ##########################################
        ax["vel"].axhspan(
            params["VEL_MIN"],
            params["VEL_MAX"],
            facecolor="gray",
            alpha=0.25,
            edgecolor=None,
        )
        if plot_size == "big":
            sc = ax["vel"].scatter(t, velocity, c=mccm, edgecolors="k", lw=scatter_lw, cmap=cm)
        else:
            sc = ax["vel"].scatter(t, velocity, c=mccm, s=8*np.ones_like(t), edgecolors="k", lw=scatter_lw, cmap=cm)
        if params["ARRAY_LABEL"] == "Hydroacoustic":
            ax["vel"].set_ylim(1.2, 1.8)
        else:
            ax["vel"].set_ylim(0.15, 0.6)
        ax["vel"].set_xlim(T1, T2)
        sc.set_clim(cax)
        if plot_size == "big":
            ax["vel"].xaxis_date()
            ax["vel"].fmt_xdata = dates.DateFormatter("%HH:%MM")
            ax["vel"].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
            ax["vel"].set_xticklabels([])
            ax["vel"].tick_params(direction="in", axis="x", top="on")
            ax["vel"].set_ylabel("Trace Velocity\n [km/s]")
        else:
            ax["vel"].set_xticks([])
            ax["vel"].set_yticks([])
        ##########################################


        ########### Plot Back-azimuths ###########
        ##########################################
        if params['AZ_MAX'] < params['AZ_MIN']:
            params['AZ_MIN'] = params['AZ_MIN']-360
            for target in params['TARGETS']:
                baz = params[target]
                if baz > 180:
                    ax["baz"].axhline(baz-360,ls='--',lw=hline_lw, color='gray',zorder=-1)
                    if plot_size == "big":
                        ax["baz"].text(t[1], baz-360, target,bbox={'facecolor':'white','edgecolor':'white','pad':0},fontsize=8,verticalalignment='center',style='italic',zorder=10)
                else:
                    ax["baz"].axhline(baz,ls='--',lw=hline_lw, color='gray',zorder=-1)
                    if plot_size == "big":
                        ax["baz"].text(t[1],baz,target,bbox={'facecolor':'white','edgecolor':'white','pad':0},fontsize=8,verticalalignment='center',style='italic',zorder=10)
            azimuth[azimuth>180]+=-360
        else:
            for target in params['TARGETS']:
                baz = params[target]
                ax["baz"].axhline(baz,ls='--',lw=hline_lw, color='gray',zorder=-1)
                if plot_size == "big":
                    ax["baz"].text(t[1],baz,target,bbox={'facecolor':'white','edgecolor':'white','pad':0},fontsize=8,verticalalignment='center',style='italic',zorder=10)
        if plot_size == "big":
            sc=ax["baz"].scatter(t,azimuth,c=mccm,edgecolors='k',lw=scatter_lw,cmap=cm,zorder=1000)
        else:
            sc=ax["baz"].scatter(t,azimuth,s=8*np.ones_like(t), c=mccm,edgecolors='k',lw=scatter_lw,cmap=cm,zorder=1000)
        ax["baz"].set_ylim(params['AZ_MIN'], params['AZ_MAX'])
        ax["baz"].set_xlim(T1,T2)
        sc.set_clim(cax)
        if plot_size == "big":
            ax["baz"].set_ylabel('Back-azimuth\n [deg]')
            ax["baz"].xaxis_date()
            ax["baz"].tick_params(axis='x',labelbottom='on')
            ax["baz"].fmt_xdata = dates.DateFormatter('%HH:%MM')
            ax["baz"].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
            ax["baz"].set_xlabel('UTC Time ['+t1.strftime('%Y-%b-%d')+']')
        else:
            ax["baz"].set_xticks([])
            ax["baz"].set_yticks([])
        ##########################################

        ########## Adjust & Save Figure ##########
        ##########################################
        if plot_size == "big":
            plt.subplots_adjust(left=0.1, right=0.9, top=0.97, bottom=0.05, hspace=0.1)
            ctop = ax["cc"].get_position().y1
            cbot = ax["baz"].get_position().y0
            cbaxes = fig.add_axes([0.91, cbot, 0.02, ctop - cbot])
            hc = plt.colorbar(sc, cax=cbaxes)
            hc.set_label(r'$M_{d}CCM$')
            filename=d2+'/'+params["ARRAY_NAME"]+'_'+t2.strftime('%Y%m%d-%H%M')+'.png'
            fig.savefig(filename,dpi=72,format='png')
            plt.close("all")
        else:
            plt.subplots_adjust(left=0,right=1,top=0.99,bottom=0.01,hspace=0.03)
            filename=d2+'/'+params["ARRAY_NAME"]+'_'+t2.strftime('%Y%m%d-%H%M')+'_thumb.png'
            fig.savefig(filename,format='png',pad_inches=0,dpi=72)
            plt.close('all')