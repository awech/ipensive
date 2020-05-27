# %% module imports
import argparse
import warnings

import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# Import the package.
import lts_array

# And the config file
import config

# Get command line arguments (if any)
parser = argparse.ArgumentParser()
# Use strftime so we always get a string out of here.
# Default if no arguments given is current time
parser.add_argument('T0', nargs='*',
                    default=[UTCDateTime.utcnow().strftime('%Y-%m-%dT%H:%M:00Z'), ])

args = parser.parse_args()
if len(args.T0) > 2:
    warnings.warn('Too many input arguments')
    parser.print_usage()
    exit(1)

ENDTIME = ''.join(args.T0)  # Time given, or current time
ENDTIME = UTCDateTime(ENDTIME)  # Convert to a UTC date/time object

# round down to the nearest 10-minute
ENDTIME = ENDTIME.replace(minute=ENDTIME.minute - (ENDTIME.minute % 10),
                          second=0,
                          microsecond=0)

STARTTIME = ENDTIME - config.duration

# %% Read in and filter data
# Array Parameters
CHAN = '*DF'  # TODO: Do these need to change, or are they all the same for all arrays?
LOC = '*'

# Filter limits
FMIN = config.f1
FMAX = config.f2

# Processing parameters
WINLEN = config.window_length
WINOVER = config.overlap

for net in config.arrays:
    NET = net['network']
    for array in net['arrays']:
        STA = array['Name']  # TODO: Or probably station identifier

        # LTS alpha parameter - subset size
        ALPHA = 0.75  # TODO: Get from array dictionary?

        print('Reading in data from IRIS')
        client = Client("IRIS")
        st = client.get_waveforms(NET, STA, LOC, CHAN,
                                  STARTTIME, ENDTIME, attach_response=True)
        st.merge(fill_value='latest')
        st.trim(STARTTIME, ENDTIME, pad='true', fill_value=0)
        st.sort()
        print(st)

        print('Removing sensitivity...')
        st.remove_sensitivity()

        stf = st.copy()
        stf.filter("bandpass", freqmin=FMIN, freqmax=FMAX, corners=2, zerophase=True)
        stf.taper(max_percentage=0.05)

        # %% Get inventory and lat/lon info
        inv = client.get_stations(network=NET, station=STA, channel=CHAN,
                                  location=LOC, starttime=STARTTIME,
                                  endtime=ENDTIME, level='channel')

        latlist = []
        lonlist = []
        staname = []
        for network in inv:
            for station in network:
                for channel in station:
                    latlist.append(channel.latitude)
                    lonlist.append(channel.longitude)
                    staname.append(channel.code)

        # Get element rijs
        rij = lts_array.getrij(latlist, lonlist)

        # Plot array coordinates as a check
        fig1 = plt.figure(1)
        plt.clf()
        plt.plot(rij[0, :], rij[1, :], 'ro')
        plt.axis('equal')
        plt.ylabel('km')
        plt.xlabel('km')
        plt.title(stf[0].stats.station)
        for i, tr in enumerate(stf):
            plt.text(rij[0, i], rij[1, i], tr.stats.location)

        # %% Run LTS array processing
        lts_vel, lts_baz, t, mdccm, stdict, sigma_tau = lts_array.ltsva(stf, rij, WINLEN, WINOVER, ALPHA)

        # %% Plotting
        # TODO: Figure out if this is the output we need, and if so, save to
        # proper directories as PNG.
        fig, axs = lts_array.lts_array_plot(stf, lts_vel, lts_baz, t, mdccm, stdict)

        d0 = config.out_web_dir + '/' + NET + '/' + STA + '/' + str(ENDTIME.year)
        d2 = d0 + '/' + '{:03d}'.format(ENDTIME.julday)
        filename = d2 + '/' + STA + '_' + ENDTIME.strftime('%Y%m%d-%H%M') + '.png'
        fig.savefig(filename, dpi = 72, format = 'png')
