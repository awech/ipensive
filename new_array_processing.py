# %% module imports
import argparse
import os
import warnings

import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.clients.fdsn import Client, header
import jinja2

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
LOC = '*'

# Filter limits
FMIN = config.f1
FMAX = config.f2

# Processing parameters
WINLEN = config.window_length
WINOVER = config.overlap

all_nets = []
all_stations = {}
for net in config.arrays:
    NET = net['network']
    NETDISP = net['display name']
    all_nets.append(NETDISP)
    all_stations[NETDISP] = []
    for array in net['arrays']:
        STA = array['id']
        STANAME = array['Name']
        all_stations[NETDISP].append(STANAME)

        CHAN = array['channel']

        # LTS alpha parameter - subset size
        ALPHA = array['Alpha']

        print('Reading in data from IRIS')
        client = Client("IRIS")
        try:

            st = client.get_waveforms(NET, STA, LOC, CHAN,
                                      STARTTIME, ENDTIME, attach_response=True)
        except header.FDSNNoDataException:
            print(f"No data retrieved for {STA}")
            continue

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

        # %% Run LTS array processing
        lts_vel, lts_baz, t, mdccm, stdict, sigma_tau = lts_array.ltsva(stf, rij, WINLEN, WINOVER, ALPHA)

        # %% Plotting
        # TODO: save to
        # proper directories as PNG.
        try:
            fig, axs = lts_array.lts_array_plot(stf, lts_vel, lts_baz, t, mdccm, stdict)
        except UnboundLocalError:
            print(f"Unable to generate plots for {STA}")
            continue

        # Generate the save path
        d2 = os.path.join(config.out_web_dir, NETDISP, STANAME, str(ENDTIME.year),
                          '{:03d}'.format(ENDTIME.julday))

        # Just for good measure, make sure it is the "real" path
        d2 = os.path.realpath(d2)

        # Make sure directory exists
        os.makedirs(d2, exist_ok = True)

        filename = os.path.join(d2, f"{STANAME}_{ENDTIME.strftime('%Y%m%d-%H%M')}.png")
        thumbnail_name = os.path.join(d2, f"{STANAME}_{ENDTIME.strftime('%Y%m%d-%H%M')}_thumb.png")

        plt.subplots_adjust(top = 0.98)

        # Adjust the colorbar positions
        # FIXME: This is ugly, but works because the two axes are always
        # added in the same order.
        # The first one is the vertical bar, the second the horizontal.
        # Would be better if we had some positive indication of which was which.
        found_first = False

        for axis in fig.axes:
            if axis not in axs:
                # This is a colorbar
                pos = axis.get_position().get_points().flatten()
                if not found_first:
                    # Vertical color bar. Move to the left.
                    found_first = True
                    pos[0] -= .03
                    pos[2] = .02
                else:
                    pass
                    # This is the horizontal color bar. Nudge it up (and to the left).
                    pos[0] -= .08
                    pos[1] += .015
                    pos[3] = .02

                axis.set_position(pos)

        fig.savefig(filename, dpi = 72, format = 'png')

        # Reconfigure plots for thumbnails
        fig.set_size_inches(4.0, 5.5)
        plt.subplots_adjust(left=0, right=0.99, bottom= 0.01, top=1.0,
                            wspace=0, hspace=0)

        for axis in fig.axes:
            if axis not in axs:
                axis.remove()  # remove colorbars
            else:
                # Remove tick marks and labels
                axis.tick_params(axis = 'both', which = 'both',
                                 bottom = False, top = False,
                                 labelbottom = False, left = False,
                                 right = False, labelleft = False)

                # Remove text labels
                for txt in axis.texts:
                    txt.remove()

        colorbar_axes = [x for x in fig.axes if x not in axs]
        for axis in colorbar_axes:
            axis.remove()

        # Lower DPI, but larger image size = smaller dots
        fig.savefig(thumbnail_name, dpi = 36, format = 'png')

# Write out the new HTML file
script_path = os.path.dirname(__file__)

with open(os.path.join(script_path, 'index.template'), 'r') as f:
    template = jinja2.Template(f.read())

html = template.render(networks = all_nets, stations = all_stations)
with open(os.path.join(config.out_web_dir, 'index.html'), 'w') as f:
    f.write(html)
