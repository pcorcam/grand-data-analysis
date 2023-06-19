# system
import os
import glob
import warnings
import argparse

# scipy
import numpy as np

# matplotlib
import matplotlib
import matplotlib.pyplot as plt

# grandlib
import grand.io.root_trees as rt

'''

##################################################
# PLOTTING SPECTROGRAMS OF TEN-SECOND GRAND DATA #
##################################################

GRANDLIB BRANCH: `dev_io_root`

This script reads out a GRAND data file in ROOT format.
It computes the spectrogram - fast fourier transform (FFT)
as a function of time and frequency - of the traces in the 
data file, where each detection unit (DU) is treated 
separately. This allows us to take a quick look at the
time dependency of the FFTs, particularly during on-site 
analyses. NOTE: this script assumes that the data file
contains ten-second (TD) data, i.e. that the time between each
event is ten seconds.

'''

##############################
# SET CUSTOM PLOT PARAMETERS #
##############################

plt.rcParams.update({'axes.labelsize': 30})
plt.rcParams.update({'xtick.labelsize': 25})
plt.rcParams.update({'ytick.labelsize': 25})
plt.rcParams.update({'axes.titlesize': 20})
plt.rcParams.update({'legend.fontsize': 20})


###########################
# DEFINE PARSER ARGUMENTS #
###########################

parser = argparse.ArgumentParser(description="Plot spectrogram of ADC traces of GRAND events\
                                              in one file.")

parser.add_argument('--path',
                    dest='path_to_data_file',
                    default='/sps/grand/data/nancay/may2023/ROOTfiles/td000400_f0005.root',
                    type=str,
                    help='Specifies the path of the ROOT data file containing the\
                          traces to plot.')

parser.add_argument("--savefig",
                    dest="savefig",
                    action='store_true',
                    help="Save plots in `../plots/` directory.")
parser.set_defaults(savefig=False)


parser.add_argument("--do_run",
                    dest="do_run",
                    action='store_true',
                    help="Run script for whole run corresponding to given data file.")
parser.set_defaults(savefig=False)

parse_args = parser.parse_args()

path_to_data_file = parse_args.path_to_data_file
savefig           = parse_args.savefig
do_run            = parse_args.do_run


###################################
# READ DATA FILE AND OBTAIN TREES #
###################################

if not os.path.exists(path_to_data_file):
    raise Exception('File not found:',path_to_data_file)

data_file = path_to_data_file.split('/')[-1]

# Nancay data explicitly has td in filename, but not the case for GP13 files
if 'td' not in data_file:
    msg = '*** Make sure this file contains ten-second data! {:} ***'.format(path_to_data_file)
    warnings.warn(msg)

# Initiate one TADC tree
df   = rt.DataFile(path_to_data_file)
tadc = df.tadc
tadc.get_entry(0)


# Check to plot for whole run or only specified file
# Adapt title and filenames for plot accordingly
if do_run:
    run = data_file.split('_')[0]
    print('Making spectrogram for whole run:',run)

    # Get all the files in the run
    data_files = sorted( glob.glob(path_to_data_file.replace(data_file,'') + run + '*.root') )
    
    title_template     = 'Run: {:}'.format(run) + ' | DU: {:} | Channel: {:}'
    plot_file_template = 'spectrogram_' + run + '_du_{:}_channel_{:}.png'

else:
    # Only use specified data file 
    data_files = [path_to_data_file]

    title_template     = 'File: {:}'.format(data_file) + ' | DU: {:} | Channel: {:}'
    plot_file_template = 'spectrogram_' + data_file.replace('.root','') + '_du_{:}_channel_{:}.png'


# Get the sampling frequency of the traces
sample_freq = tadc.adc_sampling_frequency[0] # [MHz]

# Depending on how you take data, not all ADC channels are the same.
# For GP13 data of March: X, Y, Z correspond to channels 1, 2, 3
# For Nancay data of May: X, Y, Z correspond to channels 0, 1, 2
channels = [1,2,3]
if tadc.trace_ch[0][3] == []: 
    channels = [0,1,2]


#############################################
# COMPUTE AND PLOT SPECTROGRAMS FOR EACH DU #
#############################################

du_list = tadc.get_list_of_all_used_dus()

print('File contains data from following DUs:',du_list)


# Loop over all DUs
for du in du_list:
    print('Making spectrograms for DU',du)

    fft_list_0 = []
    fft_list_1 = []
    fft_list_2 = []

    runtime_list  = []
    fft_freq_list = []

    runtime = 0 # [s]

    # Loop over all data_files
    for file in data_files:
        df   = rt.DataFile(file)
        tadc = df.tadc
        tadc.get_entry(0)

        n_entries_file = tadc.get_number_of_entries()

        msg = 'Looping over {} events in {}'.format(n_entries_file,file)
        print(msg)

        # Loop over events in data file
        for entry in range(n_entries_file):
            tadc.get_entry(entry)
            tadc.get_entry(entry)

            # Only select events with the right DU
            # NOTE: assumes that there are no coincident events
            #       i.e. only one DU (and trace etc) per event
            if du != tadc.du_id[0]:
                runtime += 10 # [s]
                continue
        
            # Get the traces
            trace_0 = tadc.trace_ch[0][channels[0]]
            trace_1 = tadc.trace_ch[0][channels[1]]
            trace_2 = tadc.trace_ch[0][channels[2]]
            
            # Get the FFTs
            fft_0 = np.abs( np.fft.rfft(trace_0) )
            fft_1 = np.abs( np.fft.rfft(trace_1) )
            fft_2 = np.abs( np.fft.rfft(trace_2) )
            
            # Add to FFT list
            fft_list_0.append(fft_0)
            fft_list_1.append(fft_1)
            fft_list_2.append(fft_2)

            # Add runtime to list and update
            runtime_list.append( np.ones(fft_0.size)*runtime/60 ) # [min]
            runtime += 10 # [s]

            # Get the frequencies of the FFT (same for all channels)
            fft_freq = np.fft.rfftfreq( len(trace_0) )*sample_freq # [MHz]

            # Add to frequency list
            fft_freq_list.append( np.ones(fft_0.size)*fft_freq )


    ########################
    # PLOT THE SPECTROGRAM #
    ########################

    ffts = [fft_list_0,fft_list_1,fft_list_2]

    # Loop over all channels
    for channel, fft in zip(channels,ffts):

        # Create figure
        fig, ax = plt.subplots(figsize=(10,8))

        # To avoid pcolormesh bug
        plt.grid(False)

        # Set colorbar scales
        # Use the mean FFT as a gauge
        vmin, vmax = 0.1*np.mean(fft), 10*np.mean(fft)

        # Make spectrogram
        cmap = ax.pcolormesh(fft_freq_list, # [MHz]
                             runtime_list, # [min]
                             fft,
                             cmap='Blues',
                             norm = matplotlib.colors.LogNorm( vmin=vmin,vmax=vmax ))
        
        # Add colorbar
        cbar = fig.colorbar(cmap, ax=ax)

        # Set axis limits
        ax.set_xlim([0,250]) # [MHz]

        # Set axis labels
        ax.set_xlabel(r'Frequency [MHz]')
        ax.set_ylabel(r'Runtime [min]')
        cbar.set_label(r'FFT [a.u.]')

        # Set figure title
        title = title_template.format(du,channel)
        ax.set_title(title)

        # Save figure
        if savefig:
            plot_dir  = '../plots/'
            plot_file = plot_file_template.format(du,channel)
            plt.savefig(plot_dir+plot_file,bbox_inches='tight',dpi=150)
            print('Saved figure at',plot_dir+plot_file)

# Show all figures
plt.show()


'''
#################
# END OF SCRIPT #
#################
'''