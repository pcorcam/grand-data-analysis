# system
import os
import glob
import argparse

# scipy
import numpy as np

# matplotlib
import matplotlib.pyplot as plt

# grandlib
import grand.io.root_trees as rt


###########################

'''

#########################################
# PLOTTING MEAN PSD OF GRAND DATA FILES #
#########################################

GRANDLIB BRANCH: `dev_io_root`

This script reads out a GRAND data file in ROOT format.
It computes the total power - the integral of the power
spectral density (PSD) over frequency - of the traces in the 
data file, where each detection unit (DU) is treated 
separately. This allows us to take a quick look at the
time dependency of the PSDs. The PSDs are computed by correctly
normalizing the fast Fourier transforms (FFTs). The total powers
are saved as npz files for each DU. NOTE: this script assumes 
that the data file contains ten-second (TD) data, i.e. that the time 
between each event for the same DU is ten seconds.

'''

##############################
# SET CUSTOM PLOT PARAMETERS #
##############################

plt.rcParams.update({'axes.labelsize': 30})
plt.rcParams.update({'xtick.labelsize': 25})
plt.rcParams.update({'ytick.labelsize': 25})
plt.rcParams.update({'axes.titlesize': 15})
plt.rcParams.update({'legend.fontsize': 20})
#plt.style.use('/pbs/home/p/pcorrea/tools/matplotlib_style_sans-serif.txt')


###########################
# DEFINE PARSER ARGUMENTS #
###########################

parser = argparse.ArgumentParser(description="Plot mean PSD of ADC traces of GRAND events\
                                              accross multiple files.")

parser.add_argument('--path',
                    dest='path_to_data_dir',
                    default='/sps/grand/data/gp13/may2023/20dB/ROOTfiles/',
                    type=str,
                    help='Specifies the path of the directory containing the\
                          ROOT files to analyse.')

parser.add_argument("--savefig",
                    dest="savefig",
                    action='store_true',
                    help="Save plots in `../plots/` directory.")
parser.set_defaults(savefig=False)

parser.add_argument("--no_show",
                    dest="show_plot",
                    action='store_true',
                    help="Show the plots.")
parser.set_defaults(no_show=False)

parser.add_argument("--file_spec",
                    dest="file_spec",
                    default='*.root',
                    type=str,
                    help="Specify which files to analyze.")

parser.add_argument('--empty_ch',
                    dest='empty_ch',
                    default=0,
                    type=int,
                    help='Specifies which of the 4 ADC channels is not branched.')


parse_args = parser.parse_args()

path_to_data_dir  = parse_args.path_to_data_dir
savefig           = parse_args.savefig
show_plot         = not parse_args.show_plot
file_spec         = parse_args.file_spec
empty_ch          = parse_args.empty_ch


#####################################
# GET ROOT FILES AND DU INFORMATION #
#####################################

if not os.path.exists(path_to_data_dir):
    raise Exception('File not found:',path_to_data_dir)

root_files = sorted( glob.glob(path_to_data_dir+file_spec) )

# Initiate TADC tree of first file to get basic info
# This assumes all files are from the same data run
# and therefore have the same run configuration
# Also ensure that the first file is not empty
number_of_dus = 0
file_idx      = 0
n_entries_tot = 0
for i, root_file in enumerate(root_files):
    tadc  = rt.TADC(root_file)
    tadc.get_entry(0)
    n_entries_tot += tadc.get_number_of_entries()

    if len( tadc.get_list_of_all_used_dus() ) > number_of_dus:
        number_of_dus = len( tadc.get_list_of_all_used_dus() )
        file_idx      = i

tadc  = rt.TADC(root_files[file_idx])
tadc.get_entry(0)

# Get the sampling frequency of the traces
sample_freq = tadc.adc_sampling_frequency[0] # [MHz]

# Compute the FFT normalization
n_samples = len( tadc.trace_ch[0][0] )
fft_freq  = np.fft.rfftfreq(n_samples) * 5e2 # [MHz]
fft_norm  = 1. / n_samples**2 / np.diff(fft_freq)[0] # [MHz^-1]

# Sometimes all four channels are enabled, but one is
# not branched
# For GP13 it is always channel 0
mask_ch = np.array(tadc.adc_enabled_channels_ch[0],dtype=bool)
if np.all(mask_ch):
    mask_ch[empty_ch] = False

# Sometimes there seems to be a bug where there are less than
# 3 enabled ADC channels
if len( mask_ch[mask_ch] ) < 3:
    mask_ch = np.ones(4,dtype=bool)
    mask_ch[empty_ch] = False

# Get some useful parameters from the data directory
path_to_data_dir_split = path_to_data_dir.split('/')
array = path_to_data_dir_split[4]
month = path_to_data_dir_split[5]
mode  = path_to_data_dir_split[6]

# Define the npz directory and file name template 
# to store the computed PSDs
npz_dir           = '/sps/grand/pcorrea/{}/{}/{}/total_power/'.format(array,month,mode)
npz_file_template = 'total_power_{}_{}_{}_du_{}.npz'.format(array,month,mode,'{}')

# Get the list of used DUs
du_list = tadc.get_list_of_all_used_dus()

msg = '\nFiles contain data from following DUs: {}'.format(du_list)
print(msg)


#########################################
# COMPUTE AND PLOT MEAN PSD FOR EACH DU #
#########################################

# Make dictionaries, that way you don't have to loop over all files
# each time you want to analyze a different DU
psd_x        = {}
psd_y        = {}
psd_z        = {}
n_entries_du = {}

for du in du_list:
    psd_x[du]        = np.zeros((n_entries_tot,fft_freq.size))
    psd_y[du]        = np.zeros((n_entries_tot,fft_freq.size))
    psd_z[du]        = np.zeros((n_entries_tot,fft_freq.size))
    n_entries_du[du] = 0

# Loop over all files in run
n_files   = len(root_files)
idx_entry = 0
for i, root_file in enumerate(root_files[:]):
    trawv = rt.TRawVoltage(root_file)

    n_entries_file = trawv.get_number_of_entries()

    msg = '\n{}/{}: Looping over {} events in {}'.format(i+1,n_files,n_entries_file,root_file)
    print(msg)

    # Loop over all entries in file
    for entry in range(n_entries_file):
        trawv.get_entry(entry)

        # Use the DU as key for the FFT dictionaries
        # NOTE: assumes that there are no coincident events
        #       i.e. only one DU (and trace etc) per event
        du = trawv.du_id[0]
    
        # Get the traces
        new_trace_ch = np.array(trawv.trace_ch[0],dtype=object)
        new_trace_ch = new_trace_ch[mask_ch]

        if len(new_trace_ch[0]) != n_samples\
        or len(new_trace_ch[1]) != n_samples\
        or len(new_trace_ch[2]) != n_samples:
            msg  = '\nWARNING for entry {}: {} samples in trace != (3,{})'.format(entry,new_trace_ch.shape,n_samples)
            msg += '\nSkipping...'
            print(msg)
            idx_entry += 1
            continue

        # Above a GPS temperature of 55ºC the ADC is unreliable
        # Mostly seen in GP13 data of March 2023
        # Sometimes this yields corrupted GPS data that is <10ºC
        # Keeping data between 10-55ºC keeps all the good data
        # Might become obsolete with furder upgrades to the prototypes
        if (trawv.gps_temp[0] > 55 or trawv.gps_temp[0] < 10)\
        and (array,month) == ('gp13','mar2023'):
            msg  = '\nWARNING for entry {}: gps_temp = {} not between 10 and 55 ºC'.format(entry,trawv.gps_temp[0])
            msg += '\nSkipping...'
            print(msg)
            idx_entry += 1
            continue

        n_entries_du[du] += 1

        trace_x = new_trace_ch[0]
        trace_y = new_trace_ch[1]
        trace_z = new_trace_ch[2]
        
        # Get the FFTs
        # Take the square for the PSD
        fft_x = np.abs( np.fft.rfft(trace_x) )**2 # [µV^2]
        fft_y = np.abs( np.fft.rfft(trace_y) )**2 # [µV^2]
        fft_z = np.abs( np.fft.rfft(trace_z) )**2 # [µV^2]
        
        # Convert to PSD and add to array contribution
        psd_x[du][idx_entry] = fft_norm * fft_x * 1e-12 # [V^2 MHz^-1]
        psd_y[du][idx_entry] = fft_norm * fft_y * 1e-12 # [V^2 MHz^-1]
        psd_z[du][idx_entry] = fft_norm * fft_z * 1e-12 # [V^2 MHz^-1]

        idx_entry += 1


for du in du_list:
    # Get rid of entries that were not filled because
    # of bad quality data, like overheating etc
    # NOTE: if there is a bad entry which is not filled,
    #       this does not take into account the time step
    #       of 10s corresponding to this gap. This means
    #       that for further analysis, each entry is assumed
    #       to be 10s after the other. In general this is not
    #       an issue.
    mask_filled = np.any(psd_x[du],axis=1)
    psd_x[du]   = psd_x[du][mask_filled]
    psd_y[du]   = psd_y[du][mask_filled]
    psd_z[du]   = psd_z[du][mask_filled]

    # Calculate the total power by integrating over frequency
    total_power_x = np.trapz(psd_x[du],x=fft_freq,axis=1) # [V^2]
    total_power_y = np.trapz(psd_y[du],x=fft_freq,axis=1) # [V^2]
    total_power_z = np.trapz(psd_z[du],x=fft_freq,axis=1) # [V^2]

    # Save total power in npz file
    npz_file = npz_file_template.format(du)
    np.savez(npz_dir+npz_file,
             du=du,
             n_entries=n_entries_du[du],
             total_power_x=total_power_x,
             total_power_y=total_power_y,
             total_power_z=total_power_z)
    msg = '\nSaved npz file: {}'.format(npz_dir+npz_file)
    print(msg)