# %% module imports
import argparse
import os
import warnings

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from obspy import UTCDateTime
from obspy.clients.fdsn import Client, header
from obspy.geodetics.base import gps2dist_azimuth
import jinja2

# Import the package.
import lts_array

# And the config file
import config


def get_volcano_backazimuth(latlist, lonlist, volcanoes):
    lon0 = np.mean(lonlist)
    lat0 = np.mean(latlist)
    for volc in volcanoes:
        if 'back_azimuth' not in volc:
            tmp = gps2dist_azimuth(lat0, lon0, volc['v_lat'], volc['v_lon'])
            volc['back_azimuth'] = tmp[1]
    return volcanoes


if __name__ == "__main__":
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
            try:
                fig, axs = lts_array.lts_array_plot(stf, lts_vel, lts_baz, t, mdccm, stdict)
            except UnboundLocalError:
                print(f"Unable to generate plots for {STA}")
                continue

            # NOTE: these are implementation dependant, and could easily change.
            backazimuth_axis = axs[2]
            velocity_axis = axs[1]

            ################ Velocity Graph ##########################
            # Tweak the y axis tick marks for the velocity plot
            v_ystart, v_yend = velocity_axis.get_ylim()
            velocity_axis.yaxis.set_ticks(np.arange(v_ystart, 0.55, 0.05))
            velocity_axis.set_ylim(top = 0.5)

            # Shade the background for the velocity area of interest
            max_vel = 0.45
            min_vel = 0.3
            velocity_axis.axhspan(min_vel, max_vel, color = "gray", zorder = -1,
                                  alpha=0.25)

            ##################### Pressure Graph ##################
            # Use thinner lines on the pressure graph
            for line in axs[0].lines:
                line.set_linewidth(0.6)

            ####################### X Axis Formatting #####################
            # Format the date/time stuff
            axs[-1].set_xlabel(f'UTC Time ({ENDTIME.strftime("%d %B %Y")})')
            axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

            ######################### Backazimuth Plot #####################
            # Add volcano azimuth lines to plots
            volcanoes = get_volcano_backazimuth(latlist, lonlist,
                                                array['volcano'])

            # decide where to put the volcano labels (horizontal position)
            # Probably overkill, I suspect we could just use start+fixed offset
            # But since I don't actually know how "long" these graphs are, I'll
            # calculate for now
            limits = backazimuth_axis.get_xlim()
            #  25 is completly arbitrary, but it seems to work nicely enough.
            label_left = limits[0] + (((limits[1] - limits[0]) / 25))

            volc_azimuth_markers = []
            for volc in volcanoes:
                # Add the line
                volc_azimuth_markers.append(backazimuth_axis.axhline(volc['back_azimuth'],
                                                                     ls = '--',
                                                                     color = "gray",
                                                                     zorder = -1))

                # And the name
                volc_azimuth_markers.append(backazimuth_axis.text(label_left,
                                                                  volc['back_azimuth'] - 6,
                                                                  volc['name'],
                                                                  bbox={'facecolor': 'white',
                                                                        'edgecolor': 'white',
                                                                        'pad': 0},
                                                                  fontsize=8,
                                                                  style='italic',
                                                                  zorder=10))

            #################### Plot Layout ##################################
            # Replace the graph title
            for txt in axs[0].texts:
                txt.remove()

            title = fig.text(.5, 0.99, f"{STANAME} Infrasound Array",
                             horizontalalignment = 'center',
                             verticalalignment = 'top')

            # Tighten up the layout
            plt.tight_layout()
            plt.subplots_adjust(top = 0.97, right = .90, bottom = 0.11)

            # Adjust the colorbar positions to not cut off
            # FIXME: This is ugly, but works because the two axes are always
            # added in the same order.
            # The first one is the vertical bar, the second the horizontal.
            # Would be better if we had some positive indication of which was which.
            vertical_colorbar = None
            horizontal_bar = None

            for axis in fig.axes:
                if axis not in axs:
                    # This is a colorbar
                    pos = axis.get_position().get_points().flatten()
                    if vertical_colorbar is None:
                        # Vertical color bar. Move to the left.
                        vertical_colorbar = axis
                        pos[0] -= .03
                        pos[1] += .01
                        pos[2] = .02
                        pos[3] -= 0.05
                    else:
                        # This is the horizontal color bar. Nudge it up (and to the left).
                        horizontal_bar = axis
                        pos[0] -= .02
                        pos[1] += .015
                        pos[2] -= 0.1
                        pos[3] = .02

                        # Make sure this one only has integer tick marks
                        _, x_max = axis.get_xlim()
                        axis.xaxis.set_ticks(np.arange(1, x_max))

                    axis.set_position(pos)
                else:
                    # Move the x axis ticks inside the plot
                    axis.tick_params(axis = 'x', direction = "in")

            ###################################################################

            # Generate the save path
            d2 = os.path.join(config.out_web_dir, NETDISP, STANAME,
                              str(ENDTIME.year),
                              '{:03d}'.format(ENDTIME.julday))

            # Just for good measure, make sure it is the "real" path.
            # Probably completly paranoid and unnecessary.
            d2 = os.path.realpath(d2)

            # Make sure directory exists
            os.makedirs(d2, exist_ok = True)

            filename = os.path.join(d2, f"{STANAME}_{ENDTIME.strftime('%Y%m%d-%H%M')}.png")
            thumbnail_name = os.path.join(d2, f"{STANAME}_{ENDTIME.strftime('%Y%m%d-%H%M')}_thumb.png")

            # Finally, save the full size image
            fig.savefig(filename, dpi = 72, format = 'png')

            # Reconfigure plots for thumbnails
            # Remove the volcano back-azimuth stuff
            for volc in volc_azimuth_markers:
                volc.remove()

            # and the plot title
            title.remove()

            # Resize down to thumbnail size and spacing
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
