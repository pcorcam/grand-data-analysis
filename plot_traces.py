# system
import os
import time
import argparse

# scipy
import numpy as np

# matplotlib
import matplotlib.pyplot as plt

# grandlib
import grand.io.root_trees as rt


'''

#######################################
# PLOTTING TRACES OF GRAND DATA FILES #
#######################################

GRANDLIB BRANCH: `dev_io_root`

This script reads out a GRAND data file in ROOT format.
It loops over the events in the file, plotting each of the ADC traces
in an interactive way. This allows us to take a quick look at traces,
particularly during on-site analyses.

'''

##############################
# SET CUSTOM PLOT PARAMETERS #
##############################

plt.rcParams.update({'axes.labelsize': 30})
plt.rcParams.update({'xtick.labelsize': 26})
plt.rcParams.update({'ytick.labelsize': 26})
plt.rcParams.update({'axes.titlesize': 10})
plt.rcParams.update({'legend.fontsize': 20})


###########################
# DEFINE PARSER ARGUMENTS #
###########################

parser = argparse.ArgumentParser(description="Plot ADC traces of GRAND events.")

parser.add_argument('--path',
                    dest='path_to_data_file',
                    default='/sps/grand/data/gp13/may2023/20dB/ROOTfiles/GRAND.TEST-RAW-10s-ChanXYZ-20dB-14dus.20230520125419.001_dat.root',
                    type=str,
                    help='Specifies the path of the ROOT data file containing the\
                          traces to plot.')

parser.add_argument('--time_sleep',
                    dest='time_sleep',
                    default=1,
                    type=float,
                    help='Specifies the sleep time between plotting traces of each entry.')

parser.add_argument('--start_entry',
                    dest='start_entry',
                    default=0,
                    type=int,
                    help='Entry of ROOT tree from where to start the plotting')

parser.add_argument('--sigma',
                    dest='sigma',
                    default=0,
                    type=int,
                    help='If sigma>0 plot lines of sigma*std.')

parser.add_argument("--savefig",
                    dest="savefig",
                    action='store_true',
                    help="Save plots in `../plots/` directory.")
parser.set_defaults(savefig=False)

parse_args = parser.parse_args()

path_to_data_file = parse_args.path_to_data_file
time_sleep        = parse_args.time_sleep
start_entry       = parse_args.start_entry
sigma             = parse_args.sigma
savefig           = parse_args.savefig


###################################
# READ DATA FILE AND OBTAIN TREES #
###################################

if not os.path.exists(path_to_data_file):
    raise Exception('File not found:',path_to_data_file)

data_file = path_to_data_file.split('/')[-1]


# Initiate TADC tree
df   = rt.DataFile(path_to_data_file)
tadc = df.tadc
#tadc = rt.TADC(path_to_data_file)


###################
# PLOT THE TRACES #
###################

# Idea is to create figure before loop, and redraw at each iteration

# Activate interactive mode
plt.ion()

# Create figure
fig, ax = plt.subplots(3,1,sharex=True,figsize=(10,8))

# Get first TADC entry
tadc.get_entry(0)

# Depending on how you take data, not all ADC channels are the same.
# For GP13 data of March: X, Y, Z correspond to channels 1, 2, 3
# For Nancay data of May: X, Y, Z correspond to channels 0, 1, 2
channels = [1,2,3]
if tadc.trace_ch[0][3] == []: 
    channels = [0,1,2]

labels = ['Channel {:}'.format(channel) for channel in channels]

# Create plt.plot instances
plt_0, = ax[0].plot(tadc.trace_ch[0][channels[0]],label=labels[0],color='b')
plt_1, = ax[1].plot(tadc.trace_ch[0][channels[1]],label=labels[1],color='m')
plt_2, = ax[2].plot(tadc.trace_ch[0][channels[2]],label=labels[2],color='r')

# Create plt.hline instances
if sigma > 0:
    label = '$\pm {}\sigma$'.format(sigma)
    hline_0_plus = ax[0].axhline(color='k',linewidth=2,linestyle='--',label=label)
    hline_1_plus = ax[1].axhline(color='k',linewidth=2,linestyle='--',label=label)
    hline_2_plus = ax[2].axhline(color='k',linewidth=2,linestyle='--',label=label)

    hline_0_minus = ax[0].axhline(color='k',linewidth=2,linestyle='--')
    hline_1_minus = ax[1].axhline(color='k',linewidth=2,linestyle='--')
    hline_2_minus = ax[2].axhline(color='k',linewidth=2,linestyle='--')

# Set axis labels
ax[2].set_xlabel('Sample number',fontsize=20)
ax[1].set_ylabel('ADC counts [LSB]',fontsize=20)

# Set figure legends
ax[0].legend(frameon=True)
ax[1].legend(frameon=True)
ax[2].legend(frameon=True)

for entry in range(tadc.get_number_of_entries())[start_entry:]:
    print(entry)

    # Load the entry in the tree
    tadc.get_entry(entry)

    # Get the traces
    trace_0 = tadc.trace_ch[0][channels[0]]
    trace_1 = tadc.trace_ch[0][channels[1]]
    trace_2 = tadc.trace_ch[0][channels[2]]

    # Get the STD of the traces
    std_0 = np.std(trace_0)
    std_1 = np.std(trace_1)
    std_2 = np.std(trace_2)

    # Plot the traces in X=NS, Y=EW, Z=UP
    plt_0.set_ydata(trace_0)
    plt_1.set_ydata(trace_1)
    plt_2.set_ydata(trace_2)

    if sigma > 0:
        hline_0_plus.set_ydata(sigma*std_0)
        hline_1_plus.set_ydata(sigma*std_1)
        hline_2_plus.set_ydata(sigma*std_2)

        hline_0_minus.set_ydata(-sigma*std_0)
        hline_1_minus.set_ydata(-sigma*std_1)
        hline_2_minus.set_ydata(-sigma*std_2)

    # Rescale y axes
    # ylim_0 = [np.min(trace_0) - 20, np.max(trace_0) + 20]
    # ylim_1 = [np.min(trace_1) - 20, np.max(trace_1) + 20]
    # ylim_2 = [np.min(trace_2) - 20, np.max(trace_2) + 20]

    if sigma < 5.5:
        lim = 7.5
    else:
        lim = sigma + 0.5

    ylim_0 = [-lim*std_0, lim*std_0]
    ylim_1 = [-lim*std_1, lim*std_1]
    ylim_2 = [-lim*std_2, lim*std_2]

    ax[0].set_ylim(ylim_0)
    ax[1].set_ylim(ylim_1)
    ax[2].set_ylim(ylim_2)

    # Set figure title
    title = r'File: {:} | Entry: {:} | Event index: {:} | DU: {:}'
    title = title.format(data_file,entry,tadc.event_id[0],tadc.du_id[0])
    ax[0].set_title(title)

    # Draw figure
    fig.canvas.draw()

    # Save figure if specified
    if savefig:
        plot_dir  = '/pbs/home/p/pcorrea/GRAND/plots/traces/'
        extension = '_entry_{}.png'.format(entry)
        plot_file = data_file.replace('.root',extension)
        plt.savefig(plot_dir+plot_file,bbox_inches='tight',dpi=150)
        print('Saved figure at',plot_dir+plot_file)

    # Flush events for next iterations
    fig.canvas.flush_events()
    
    # Sleep for some time to show the figure
    time.sleep(time_sleep)

'''
#################
# END OF SCRIPT #
#################
'''