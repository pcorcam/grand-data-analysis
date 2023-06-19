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

###################################
# GET THE RMS OF GRAND DATA FILES #
###################################

GRANDLIB BRANCH: `dev_io_root`

This script reads out a GRAND data file in ROOT format.
It computes the root mean square (RMS) of the traces
in the data file, where each detection unit (DU) is treated
separately. The RMS is computed as function of time, and 
saved in an npz file for each DU.

'''


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
n_samples = len( tadc.trace_ch[0][0] )


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
npz_dir           = '/sps/grand/pcorrea/{}/{}/{}/rms/'.format(array,month,mode)
npz_file_template = 'rms_{}_{}_{}_du_{}.npz'.format(array,month,mode,'{}')

# Get the list of used DUs
du_list = tadc.get_list_of_all_used_dus()

msg = '\nFiles contain data from following DUs: {}'.format(du_list)
print(msg)


###############################################
# COMPUTE RMS AS FUNCTION OF TIME FOR EACH DU #
###############################################

# Make dictionaries, that way you don't have to loop over all files
# each time you want to analyze a different DU
rms_x    = {}
rms_y    = {}
rms_z    = {}
gps_time = {}

for du in du_list:
    rms_x[du]    = np.zeros((n_entries_tot))
    rms_y[du]    = np.zeros((n_entries_tot))
    rms_z[du]    = np.zeros((n_entries_tot))
    gps_time[du] = np.zeros((n_entries_tot))

# Loop over all files in run
n_files = len(root_files)
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
            continue
            
        trace_x = new_trace_ch[0]
        trace_y = new_trace_ch[1]
        trace_z = new_trace_ch[2]
        
        # Compute the RMS
        rms_x[du][idx_entry] = np.sqrt( np.mean( trace_x**2 ) )
        rms_y[du][idx_entry] = np.sqrt( np.mean( trace_y**2 ) )
        rms_z[du][idx_entry] = np.sqrt( np.mean( trace_z**2 ) )

        # Get the GPS timestamp
        gps_time[du][idx_entry] = trawv.gps_time[0]

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
    mask_filled  = np.where(rms_x[du]>0)
    rms_x[du]    = rms_x[du][mask_filled]
    rms_y[du]    = rms_y[du][mask_filled]
    rms_z[du]    = rms_z[du][mask_filled]
    gps_time[du] = gps_time[du][mask_filled]

    # Save total power in npz file
    npz_file = npz_file_template.format(du)
    np.savez(npz_dir+npz_file,
             du=du,
             rms_x=rms_x[du],
             rms_y=rms_y[du],
             rms_z=rms_z[du],
             gps_time=gps_time[du])
    msg = '\nSaved npz file: {}'.format(npz_dir+npz_file)
    print(msg)


