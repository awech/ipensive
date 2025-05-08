import os
import sys
import numpy as np
import pandas as pd
from obspy import Stream, UTCDateTime, read_inventory
from obspy.clients.earthworm import Client
from obspy.clients.fdsn import Client as FDSNClient
from obspy.geodetics.base import gps2dist_azimuth
import yaml

####### plotting imports #######
import matplotlib as m
m.use('Agg')
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib import dates
from copy import deepcopy
from collections import Counter
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

def write_ascii_file(t2, t, pressure, azimuth, velocity, mccm, rms, name, config):


    t1=t2-config[name]["DURATION"]
    tmp_name=name.replace(' ','_')

    d0=config["OUT_ASCII_DIR"]+'/'+tmp_name
    if not os.path.exists(d0):
        os.mkdir(d0)

    subfolder=d0+'/{}'.format(t1.strftime('%Y-%m'))
    if not os.path.exists(subfolder):
        os.mkdir(subfolder)

    
    filename=subfolder+'/'+tmp_name+'_'+t1.strftime('%Y-%m-%d')+'.txt'

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


def write_valve_file(t2, t, pressure, azimuth, velocity, mccm, rms, name, config):

    A=pd.DataFrame({'TIMESTAMP':t,
                 'CHANNEL':name,
                 'Azimuth':azimuth,
                 'Velocity':1000*velocity,
                 'MCCM':mccm,
                 'Pressure':pressure,
                 'rms':rms})

    # A=A[['TIMESTAMP','CHANNEL','Azimuth','Velocity','MCCM','Pressure','rms']]

    # A['Velocity']=1000*A['Velocity']

    filename=config["OUT_VALVE_DIR"]+'/'+name+'_'+t2.strftime('%Y%m%d-%H%M')+'.txt'

    A.to_csv(filename,index=False,header=True,sep=',',float_format='%.3f')


def plot_results(t1, t2, t, st, mccm, velocity, azimuth, lts_dict, config, params):

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

    if params["PLOT_MCCM"]:
        ax_list = [["wave"], ["cc"], ["vel"], ["baz"], ["stas"]]
        hr_list = [1, 1, 1, 1, 0.66]
    else:
        ax_list = [["wave"], ["vel"], ["baz"], ["stas"]]
        hr_list = [1, 1, 1, 0.66]
    size = (8, 10.5)
    trace_lw = 0.6
    scatter_lw = 0.3
    s_dot = 36
    hline_lw = 1
    wm_font = 18
    for plot_size in ["big", "small"]:
        if plot_size == "small":
            size = (2.1, 2.75)
            trace_lw = 0.1
            scatter_lw = 0.1
            hline_lw = 0.25
            s_dot = 8
            wm_font = 8

        fig, ax = plt.subplot_mosaic(
            ax_list,
            figsize=size,
            height_ratios=hr_list
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
        if params["PLOT_MCCM"]:
            sc = ax["cc"].scatter(t, mccm, c=mccm, s=s_dot, edgecolors="k", lw=scatter_lw, cmap=cm)
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

        sc = ax["vel"].scatter(t, velocity, c=mccm, s=s_dot, edgecolors="k", lw=scatter_lw, cmap=cm)
        
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
        box_style = {'facecolor':'white','edgecolor':'white','pad':0}
        az_min = deepcopy(params['AZ_MIN'])
        az_max = deepcopy(params['AZ_MAX'])
        tmp_azimuth = deepcopy(azimuth)
        if params['AZ_MAX'] < params['AZ_MIN']:
            az_min = az_min-360
            for target in params['TARGETS']:
                baz = params[target]
                if baz > 180:
                    ax["baz"].axhline(baz-360, ls='--', lw=hline_lw, color='gray', zorder=-1)
                    if plot_size == "big":
                        ax["baz"].text(t[1], baz-360, target, bbox=box_style, fontsize=8, va='center', style='italic', zorder=10)
                else:
                    ax["baz"].axhline(baz, ls='--', lw=hline_lw, color='gray', zorder=-1)
                    if plot_size == "big":
                        ax["baz"].text(t[1], baz, target, bbox=box_style, fontsize=8, va='center', style='italic', zorder=10)
            tmp_azimuth[tmp_azimuth>180]+=-360
        else:
            for target in params['TARGETS']:
                baz = params[target]
                ax["baz"].axhline(baz, ls='--', lw=hline_lw, color='gray', zorder=-1)
                if plot_size == "big":
                    ax["baz"].text(t[1], baz, target, bbox=box_style, fontsize=8, va='center', style='italic', zorder=10)

        sc=ax["baz"].scatter(t, tmp_azimuth, c=mccm, s=s_dot, edgecolors='k', lw=scatter_lw, cmap=cm, zorder=1000)
        ax["baz"].set_ylim(az_min, az_max)
        ax["baz"].set_xlim(T1,T2)
        sc.set_clim(cax)
        if plot_size == "big":
            ax["baz"].xaxis_date()
            ax["baz"].fmt_xdata = dates.DateFormatter("%HH:%MM")
            ax["baz"].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
            ax["baz"].set_xticklabels([])
            ax["baz"].tick_params(direction="in", axis="x", top="on")
            ax["baz"].set_ylabel("Trace Velocity\n [km/s]")
        else:
            ax["baz"].set_xticks([])
            ax["baz"].set_yticks([])

        ########## Plot Dropped Channels #########
        ##########################################
        if len(lts_dict) == 0:
            txt_str = "Not enough channels for LTS"
            ax["stas"].text(0.5, 0.5, txt_str, transform=ax["stas"].transAxes, color='grey', alpha=0.7, fontsize=wm_font, va="center", ha="center")
            n = len(st)
        else:
            ndict = deepcopy(lts_dict)
            n = ndict['size']
            ndict.pop('size', None)
            tstamps = list(ndict.keys())
            tstampsfloat = [float(ii) for ii in tstamps]
            cm2 = plt.get_cmap('binary', (n-1))
            for jj in range(len(tstamps)):
                z = Counter(list(ndict[tstamps[jj]]))
                keys, vals = z.keys(), z.values()
                keys, vals = np.array(list(keys)), np.array(list(vals))
                pts = np.tile(tstampsfloat[jj], len(keys))
                sc_stas = ax["stas"].scatter(
                    pts,
                    keys,
                    c=vals,
                    s=s_dot,
                    edgecolors="k",
                    lw=scatter_lw,
                    cmap=cm2,
                    vmin=0.5,
                    vmax=n - 0.5,
                )
        ax["stas"].set_ylim(0, n+1)
        ax["stas"].invert_yaxis()
        ax["stas"].set_yticks(np.arange(1, n+1))
        ax["stas"].set_xlim(T1,T2)
        if plot_size == "big":
            ax["stas"].set_ylabel('Element [#]')
            ax["stas"].xaxis_date()
            ax["stas"].tick_params(axis='x',labelbottom='on')
            ax["stas"].fmt_xdata = dates.DateFormatter('%HH:%MM')
            ax["stas"].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
            ax["stas"].set_xlabel('UTC Time ['+t1.strftime('%Y-%b-%d')+']')
        else:
            ax["stas"].set_xticks([])
            ax["stas"].set_yticks([])
        ##########################################

        ########## Adjust & Save Figure ##########
        ##########################################
        if plot_size == "big":
            plt.subplots_adjust(left=0.1, right=0.9, top=0.97, bottom=0.05, hspace=0.1)

            ax_str = "cc" if params["PLOT_MCCM"] else "vel"
            ctop = ax[ax_str].get_position().y1
            cbot = ax["baz"].get_position().y0
            cbaxes_mccm = fig.add_axes([0.91, cbot, 0.02, ctop - cbot])
            hc = plt.colorbar(sc, cax=cbaxes_mccm)
            hc.set_label(r'$M_{d}CCM$')
            if n > 3 and params["LTS_ALPHA"] < 1:
                ctop = ax["stas"].get_position().y1
                cbot = ax["stas"].get_position().y0
                cbaxes_stas = fig.add_axes([0.91, cbot, 0.02, ctop - cbot])
                hc_stas = plt.colorbar(sc_stas, cax=cbaxes_stas)
                hc_stas.set_label(r"# Dropped Pairs")
                hc_stas.set_ticks(np.arange(1, n))
            
            filename=d2+'/'+params["ARRAY_NAME"]+'_'+t2.strftime('%Y%m%d-%H%M')+'.png'
            fig.savefig(filename,dpi=72,format='png')
            plt.close("all")
        else:
            plt.subplots_adjust(left=0,right=1,top=0.99,bottom=0.01,hspace=0.03)
            filename=d2+'/'+params["ARRAY_NAME"]+'_'+t2.strftime('%Y%m%d-%H%M')+'_thumb.png'
            fig.savefig(filename,format='png',pad_inches=0,dpi=72)
            plt.close('all')
