import logging
from pathlib import Path
from copy import deepcopy
from collections import Counter
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
from matplotlib import dates
from matplotlib import rcParams

rcParams.update({"font.size": 10})  # Set default font size for plots
m.use("Agg")  # Use a non-interactive backend for matplotlib
################################

my_log = logging.getLogger(__name__)


def missing_elements_adjust(st, skip_chans, array):
    """
    Adjust the indices of an array to account for missing channels.

    Args:
        st (Stream): ObsPy Stream object containing traces.
        skip_chans (list): List of channels to skip.
        array (numpy.ndarray): Array of indices to adjust.

    Returns:
        tuple: Adjusted array and a list of missing indices.
    """
    skip_chans.sort()  # Ensure the skipped channels are sorted
    missing_inds = []  # List to store indices of missing channels
    for chan in skip_chans:
        # Find the index of the missing channel in the stream
        ind = np.where(np.array([tr.id for tr in st]) == chan)[0]
        # Increment indices in the array that are greater than the missing index
        array[array > ind] += 1
        # Add the missing index (adjusted by +1) to the list
        missing_inds.append(ind + 1)
    return array, missing_inds


def update_axes_and_ticks(ax, ylabel):
    """
    Update the x-axis of a matplotlib Axes object to display date-formatted ticks and set the y-axis label.

    Parameters:
        ax (matplotlib.axes.Axes): The axes object to update.
        ylabel (str): The label to set for the y-axis.

    Effects:
        - Formats the x-axis to display dates in "%H:%M" format.
        - Removes x-axis tick labels.
        - Sets the direction of x-axis ticks to "in" and enables ticks on the top of the axis.
        - Sets the y-axis label to the provided string.
    """
    ax.xaxis_date()
    ax.fmt_xdata = dates.DateFormatter("%HH:%MM")
    ax.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
    ax.set_xticklabels([])
    ax.tick_params(direction="in", axis="x", top="on")
    ax.set_ylabel(ylabel)


def plot_waveform(ax, T1, T2, st, array_params, plot_params):
    """
    Plot the waveform data for the given stream.

    Args:
        ax (matplotlib.axes.Axes): The axes to plot on.
        T1 (matplotlib datenum): The start value for the x-axis limit.
        T2 (matplotlib datenum): The end value for the x-axis limit.
        st (Stream): ObsPy Stream object containing traces.
        array_params (dict): Parameters for the array being processed.
        plot_params (dict): Parameters for the plot, including size and line widths.
    """

    # Generate a time vector for plotting waveforms
    tvec = np.linspace(
        dates.date2num(st[0].stats.starttime.datetime),
        dates.date2num(st[0].stats.endtime.datetime),
        len(st[0].data),
    )
    
    ax.set_title(array_params["ARRAY_NAME"] + " " + array_params["ARRAY_LABEL"] + " Array")
    ax.plot(tvec, st[0].data, "k", linewidth=plot_params["trace_lw"])
    ax.axis("tight")
    ax.set_xlim(T1, T2)
    ymax = np.abs(list(ax.get_ylim())).max()
    ax.set_ylim(-ymax, ymax)
    if plot_params["plot_size"] == "big":
        update_axes_and_ticks(ax, "Pressure [Pa]")
        ax2 = ax.twinx()
        ax2.set_yticks([])
        ax2.set_ylabel(
            f"{array_params['FREQMIN']:.1f} - {array_params['FREQMAX']:.1f} Hz",
            labelpad=6,
        )
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    ##########################################


def plot_cc_values(ax, T1, T2, t, mccm, array_params, plot_params):
    """
    Plots MCCM (Multi-Channel Cross-Correlation Matrix) values on the given matplotlib axis with color mapping and formatting.

    Parameters:
        ax (matplotlib.axes.Axes): The axis on which to plot.
        T1 (matplotlib datenum): The start time for the x-axis limit.
        T2 (matplotlib datenum): The end time for the x-axis limit.
        t (array-like): Time (datenum) values.
        mccm (array-like): The MCCM values to plot on the y-axis and for color mapping.
        array_params (dict): Dictionary containing plot parameters, must include "MCTHRESH" for threshold line.
        plot_params (dict): Dictionary containing plot parameters, including size and line widths.

    Returns:
        matplotlib.collections.PathCollection: The scatter plot object for further customization or colorbar attachment.
    """

    # Scatter plot of MCCM values with color mapping
    sc = ax.scatter(t, mccm, c=mccm, s=plot_params["s_dot"], edgecolors="k", lw=plot_params["scatter_lw"], cmap=plot_params["cm"])
    # Add a horizontal line for the MCCM threshold
    ax.axhline(array_params["MCTHRESH"], ls="--", lw=plot_params["hline_lw"], color="gray")
    # Adjust axis limits and scaling
    ax.axis("tight")
    ax.set_xlim(T1, T2)
    ax.set_ylim(plot_params["cax"])
    sc.set_clim(plot_params["cax"])
    # Configure axis labels and formatting for larger plots
    if plot_params["plot_size"] == "big":
        update_axes_and_ticks(ax, r"$M_{d}CCM$")
    else:
        # Remove ticks and labels for smaller plots
        ax.set_xticks([])
        ax.set_yticks([])

    return sc


def plot_trace_velocities(ax, T1, T2, t, velocity, mccm, array_params, plot_params):
    """
    Plots trace velocities on the given matplotlib axis with color mapping and formatting.

    Parameters:
        ax (matplotlib.axes.Axes): The axis on which to plot.
        T1 (matplotlib datenum): The start time for the x-axis limit.
        T2 (matplotlib datenum): The end time for the x-axis limit.
        t (array-like): Time (datenum) values.
        velocity (array-like): The trace velocities to plot.
        mccm (array-like): The MCCM values for color mapping.
        array_params (dict): Dictionary containing plot parameters.
        plot_params (dict): Dictionary containing plot parameters, including size and line widths.

    Returns:
        matplotlib.collections.PathCollection: The scatter plot object for further customization or colorbar attachment.
    """
    
    # Highlight the velocity range of interest with a shaded region
    ax.axhspan(
        array_params["VEL_MIN"],
        array_params["VEL_MAX"],
        facecolor="gray",
        alpha=0.25,
        edgecolor=None,
        )

    # Scatter plot of trace velocities with color mapping based on MCCM values
    sc = ax.scatter(t, velocity, c=mccm, s=plot_params["s_dot"], edgecolors="k", lw=plot_params["scatter_lw"], cmap=plot_params["cm"])
    # Set y-axis limits based on the array type
    if array_params["ARRAY_LABEL"] == "Hydroacoustic":
        ax.set_ylim(1.2, 1.8)  # Typical range for hydroacoustic arrays
    else:
        ax.set_ylim(0.15, 0.6)  # Typical range for other arrays
    # Set x-axis limits to match the time range
    ax.set_xlim(T1, T2)
    # Set the color limits for the scatter plot
    sc.set_clim(plot_params["cax"])
    # Configure axis labels and formatting for larger plots
    if plot_params["plot_size"] == "big":
        update_axes_and_ticks(ax, "Trace Velocity\n [km/s]")
    else:
        # Remove ticks and labels for smaller plots
        ax.set_xticks([])
        ax.set_yticks([])
        ##########################################
    return sc


def plot_back_azimuths(ax, T1, T2, t, azimuth, mccm, array_params, plot_params):
    """
    Plot back-azimuth values on the given matplotlib axis with color mapping and formatting.

    Args:
        ax (dict): Dictionary of matplotlib axes, must contain "baz".
        T1 (matplotlib datenum): Start time for the x-axis limit.
        T2 (matplotlib datenum): End time for the x-axis limit.
        t (array-like): Time (datenum) values.
        azimuth (array-like): Back-azimuth values to plot.
        mccm (array-like): MCCM values for color mapping.
        array_params (dict): Dictionary containing plot parameters, must include 'AZ_MIN', 'AZ_MAX', 'TARGETS'.
        plot_params (dict): Dictionary containing plot parameters, including size and line widths.

    Returns:
        matplotlib.collections.PathCollection: The scatter plot object for further customization or colorbar attachment.
    """

    
    # Configure the style for the text boxes
    box_style = {'facecolor': 'white', 'edgecolor': 'white', 'pad': 0}
    # Deepcopy azimuth min/max to avoid modifying input
    az_min = deepcopy(array_params['AZ_MIN'])
    az_max = deepcopy(array_params['AZ_MAX'])
    tmp_azimuth = deepcopy(azimuth)

    # If AZ_MAX < AZ_MIN, adjust azimuth range and values for plotting
    if array_params['AZ_MAX'] < array_params['AZ_MIN']:
        az_min = az_min - 360
        # Loop through targets and adjust back-azimuths
        for target in array_params['TARGETS']:
            # If target is a dict, extract its key
            if isinstance(target, dict):
                target = list(target.keys())[0]
            baz = array_params[target]
            if baz > 180:
                # Plot horizontal lines and labels for targets above 180
                ax.axhline(baz - 360, ls='--', lw=plot_params["hline_lw"], color='gray', zorder=-1)
                if plot_params["plot_size"] == "big":
                    ax.text(t[1], baz - 360, target, bbox=box_style, fontsize=8, va='center', style='italic', zorder=10)
            else:
                # Plot horizontal lines and labels for targets below 180
                ax.axhline(baz, ls='--', lw=plot_params["hline_lw"], color='gray', zorder=-1)
                if plot_params["plot_size"] == "big":
                    ax.text(t[1], baz, target, bbox=box_style, fontsize=8, va='center', style='italic', zorder=10)
        # Adjust azimuth values above 180 for plotting
        tmp_azimuth[tmp_azimuth > 180] += -360
    else:
        # Plot horizontal lines and labels for all targets
        for target in array_params['TARGETS']:
            if type(target) is dict:
                target = list(target.keys())[0]
            baz = array_params[target]
            ax.axhline(baz, ls='--', lw=plot_params["hline_lw"], color='gray', zorder=-1)
            if plot_params["plot_size"] == "big":
                ax.text(t[1], baz, target, bbox=box_style, fontsize=8, va='center', style='italic', zorder=10)

    # Scatter plot for back-azimuth values, colored by MCCM
    sc = ax.scatter(t, tmp_azimuth, c=mccm, s=plot_params["s_dot"], edgecolors='k', lw=plot_params["scatter_lw"], cmap=plot_params["cm"], zorder=1000)
    ax.set_ylim(az_min, az_max)
    ax.set_xlim(T1, T2)
    sc.set_clim(plot_params["cax"])

    if plot_params["plot_size"] == "big":
        # Configure x-axis for larger plots
        update_axes_and_ticks(ax, "Back-Azimuth\n [deg]")
    else:
        # Remove ticks for smaller plots
        ax.set_xticks([])
        ax.set_yticks([])

    return sc


def plot_lts_dropped_channels(ax, T1, T2, t, st, lts_dict, skip_chans, plot_params):
    """
    Plot dropped channels (LTS outliers) on the provided axis.

    Args:
        ax (dict): Dictionary of matplotlib axes, must contain "stas".
        T1 (matplotlib datenum): The start time for the x-axis limit.
        T2 (matplotlib datenum): The end time for the x-axis limit.
        t (numpy.ndarray): Array of timestamps for the results.
        st (Stream): ObsPy Stream object containing traces.
        lts_dict (dict): Dictionary of LTS (Least Trimmed Squares) results.
        skip_chans (list): List of channels to skip.
        plot_params (dict): Dictionary containing plot parameters, including size and line widths.

    Returns:
        sc_stas: scatter plot object for dropped channels (or None if not plotted)
    """
    
    sc_stas = None
    missing_inds = []
    ndict = {}
    n = len(st)

    if len(lts_dict) == 0:
        # Handle case where there are not enough channels for LTS
        txt_str = "Not enough channels for LTS"
        ax.text(
            0.5, 0.5, txt_str,
            transform=ax.transAxes,
            color='grey', alpha=0.7,
            fontsize=plot_params["wm_font"], va="center", ha="center"
        )
        # Identify missing channels
        for chan in skip_chans:
            ind = np.where(np.array([tr.id for tr in st]) == chan)[0] + 1
            missing_inds.append(ind)
    else:
        # Process LTS dictionary to plot dropped channels
        ndict = deepcopy(lts_dict)
        n = ndict['size']
        ndict.pop('size', None)
        tstamps = list(ndict.keys())
        tstampsfloat = [float(ii) for ii in tstamps]
        cm2 = plt.get_cmap('binary', (n - 1))
        for jj in range(len(tstamps)):
            z = Counter(list(ndict[tstamps[jj]]))
            keys, vals = z.keys(), z.values()
            keys, vals = np.array(list(keys)), np.array(list(vals))
            pts = np.tile(tstampsfloat[jj], len(keys))
            keys, missing_inds = missing_elements_adjust(st, skip_chans, keys)
            # Scatter plot for dropped channels
            sc_stas = ax.scatter(
                pts,
                keys,
                c=vals,
                s=plot_params["s_dot"],
                edgecolors="k",
                lw=plot_params["scatter_lw"],
                cmap=cm2,
                vmin=0.5,
                vmax=n - 0.5,
            )
    # Configure y-axis for station labels
    ax.set_ylim(0, len(st) + 1)
    ax.set_yticks(np.arange(1, len(st) + 1))
    ax.set_yticklabels([f"{tr.stats.station}.{tr.stats.location}" for tr in st], fontsize=8)
    ax.set_xlim(T1, T2)

    if plot_params["plot_size"] == "big":
        # Add scatter plot for missing channels
        if n > 3 and ndict:
            for m_ind in missing_inds:
                ax.scatter(t, m_ind * np.ones(t.shape[0]), marker=".", color="indianred", s=10)

        # Configure x-axis for larger plots
        update_axes_and_ticks(ax, "")

    else:
        # Remove ticks for smaller plots
        ax.set_xticks([])
        ax.set_yticks([])
    ax.invert_yaxis()
    ##########################################

    return sc_stas


def add_lts_colorbar(ax, fig, sc_stas, lts_dict, plot_params):
    """Add a colorbar for the LTS (Least Trimmed Squares) results.

    Args:
        ax (matplotlib.axes.Axes): The axes to add the colorbar to.
        fig (matplotlib.figure.Figure): The figure containing the axes.
        sc_stas (matplotlib.collections.PathCollection): The scatter plot object for dropped channels.
        lts_dict (dict): Dictionary of LTS results.
        plot_params (dict): Dictionary containing plot parameters.
    """

    if sc_stas is not None:
        ctop = ax.get_position().y1
        cbot = ax.get_position().y0
        cbaxes_stas = fig.add_axes([0.91, cbot, 0.02, ctop - cbot])
        hc_stas = plt.colorbar(sc_stas, cax=cbaxes_stas)
        x = plot_params["lts_alpha"]
        hc_stas.set_label("# Dropped Pairs\n"+rf"$\alpha={x}$")
        hc_stas.set_ticks(np.arange(1, lts_dict['size']))
    elif lts_dict:
        txt_str = "No dropped channels"
        ax.text(0.5, 0.5, txt_str, transform=ax.transAxes, color='grey', alpha=0.7, fontsize=plot_params["wm_font"], va="center", ha="center")


def add_mccm_colorbar(ax1, ax2, fig, sc):
    """
    Add a colorbar for the MCCM (Multi-Channel Cross-Matching) results.

    Args:
        ax1 (matplotlib.axes.Axes): The axes for the first subplot.
        ax2 (matplotlib.axes.Axes): The axes for the second subplot.
        fig (matplotlib.figure.Figure): The figure containing the axes.
        sc (matplotlib.collections.PathCollection): The scatter plot object for MCCM results.
    """
    
    ctop = ax1.get_position().y1
    cbot = ax2.get_position().y0
    cbaxes_mccm = fig.add_axes([0.91, cbot, 0.02, ctop - cbot])
    hc = plt.colorbar(sc, cax=cbaxes_mccm)
    hc.set_label(r'$M_{d}CCM$')


def save_figure(fig, config, array_params, t2, plot_size):
    """
    Save the figure to a file.

    Args:
        fig (matplotlib.figure.Figure): The figure to save.
        config (dict): Configuration dictionary.
        array_params (dict): Array parameters.
        t2 (UTCDateTime): End time of the data window.
        plot_size (str): Size of the plot ("big" or "small").
    """
    
    # Define output directories for plots
    out_dir = Path(config["OUT_WEB_DIR"]) / array_params["NETWORK_NAME"] / array_params["ARRAY_NAME"] / str(t2.year)
    out_dir = out_dir / '{:03d}'.format(t2.julday)

    if plot_size == "big":
        filename = out_dir / (array_params["ARRAY_NAME"] + '_' + t2.strftime('%Y%m%d-%H%M') + '.png')
        fig.savefig(filename, dpi=72, format='png')

    elif plot_size == "small":
        filename = out_dir / (array_params["ARRAY_NAME"] + '_' + t2.strftime('%Y%m%d-%H%M') + '_thumb.png')
        fig.savefig(filename, format='png', pad_inches=0, dpi=72)

    plt.close('all')


def plot_results(t1, t2, t, st, mccm, velocity, azimuth, lts_dict, skip_chans, config, array_params, plot_size):
    """
    Generate and save plots for the results of array processing.

    Args:
        t1 (UTCDateTime): Start time of the data window.
        t2 (UTCDateTime): End time of the data window.
        t (numpy.ndarray): Array of timestamps for the results.
        st (Stream): ObsPy Stream object containing traces.
        mccm (numpy.ndarray): MCCM values for the results.
        velocity (numpy.ndarray): Trace velocities for the results.
        azimuth (numpy.ndarray): Back-azimuth values for the results.
        lts_dict (dict): Dictionary of LTS (Least Trimmed Squares) results.
        config (dict): Configuration dictionary.
        array_params (dict): Parameters for the array being processed.
        skip_chans (list): List of channels to skip.

    Returns:
        None
    """

    my_log.info(f"Making {plot_size} plot...")

    T1 = dates.date2num(t1.datetime)
    T2 = dates.date2num(t2.datetime)

    plot_params = {"plot_size": plot_size}
    # Define colormap and color axis limits for MCCM values
    plot_params["cm"] = "RdYlBu_r"
    plot_params["cax"] = 0.2, 1  # Colorbar range for MCCM values
    plot_params["lts_alpha"] = array_params["LTS_ALPHA"]  # Alpha value for LTS processing

    # Set default plot size and styling parameters
    if plot_size == "big":
        size = (8, 10.5)
        plot_params["trace_lw"] = 0.6
        plot_params["scatter_lw"] = 0.3
        plot_params["s_dot"] = 36
        plot_params["hline_lw"] = 1
        plot_params["wm_font"] = 18
        left, right, top, bottom, hspace = 0.1, 0.9, 0.97, 0.05, 0.1
    elif plot_size == "small":
        size = (2.1, 2.75)
        plot_params["trace_lw"] = 0.1
        plot_params["scatter_lw"] = 0.1
        plot_params["hline_lw"] = 0.4
        plot_params["s_dot"] = 8
        plot_params["wm_font"] = 8
        left, right, top, bottom, hspace = 0, 1, 0.99, 0.01, 0.03

    # Determine the layout of the plot based on whether MCCM is being plotted
    if array_params["PLOT_MCCM"]:
        ax_list = [["wave"], ["cc"], ["vel"], ["baz"], ["stas"]]
        hr_list = [1, 1, 1, 1, 0.66]
    else:
        ax_list = [["wave"], ["vel"], ["baz"], ["stas"]]
        hr_list = [1, 1, 1, 0.66]

    # Create the figure and axes using a mosaic layout
    fig, ax = plt.subplot_mosaic(
        ax_list,
        figsize=size,
        height_ratios=hr_list,
        gridspec_kw={"left": left, "right": right, "top": top, "bottom": bottom, "hspace": hspace}
    )

    ############ Plot waveforms ##############
    ##########################################
    plot_waveform(ax["wave"], T1, T2, st, array_params, plot_params)

    ############ Plot cc values ##############
    ##########################################
    # if MCCM (Mean Cross-Correlation Metric) plotting is enabled
    if array_params["PLOT_MCCM"]:
        sc = plot_cc_values(
            ax["cc"], T1, T2, t, mccm, array_params, plot_params
        )

    ########## Plot Trace Velocities #########
    ##########################################
    sc = plot_trace_velocities(
        ax["vel"], T1, T2, t, velocity, mccm, array_params, plot_params
    )

    ########### Plot Back-azimuths ###########
    ##########################################
    sc = plot_back_azimuths(
        ax["baz"], T1, T2, t, azimuth, mccm, array_params, plot_params
    )

    ########## Plot Dropped Channels #########
    ##########################################
    sc_stas = plot_lts_dropped_channels(
        ax["stas"], T1, T2, t, st, lts_dict, skip_chans, plot_params
    )

    ############## Add Colorbars #############
    ##########################################
    if plot_size == "big":
        add_lts_colorbar(ax["stas"], fig, sc_stas, lts_dict, plot_params)
        ax_str = "cc" if array_params["PLOT_MCCM"] else "vel"
        add_mccm_colorbar(ax[ax_str], ax["baz"], fig, sc)

        ax_str = ax_list[-1][0]
        ax[ax_str].xaxis_date()
        ax[ax_str].tick_params(axis='x', labelbottom='on', )
        ax[ax_str].fmt_xdata = dates.DateFormatter('%HH:%MM')
        ax[ax_str].xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
        ax[ax_str].set_xlabel('UTC Time [' + t1.strftime('%Y-%b-%d') + ']')

    save_figure(fig, config, array_params, t2, plot_size)
