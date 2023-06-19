# system
import os
import glob
import argparse
import datetime

# scipy
import numpy as np

# grandlib
import grand.io.root_trees as rt


###########################

'''

################################################
# SEARCHING FOR TRANSIENTS IN GRAND DATA FILES #
################################################

GRANDLIB BRANCH: `dev_io_root`

This script reads out a GRAND data file in ROOT format.
It searches for transient pulses in the time traces of
the data file, where each detection unit (DU) is treated
separately. Here, a transient is defined as a pulse that
exceeds 5 times the standard deviation (STD) of the trace,
called 5sigma from this point onwards. For each trace, 
the number of 5sigma transients is computed, as well as 
the number of times this 5sigma threshold is exceeded
during the transient. The maximum peak position of the 
transient is also stored. All information is saved in 
an npz file for each DU.

'''


###########################
# DEFINE PARSER ARGUMENTS #
###########################

parser = argparse.ArgumentParser(description="Search for transient in data traces and\
                                              save information in npz files per DU.")

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

parser.add_argument('--thresh',
                    dest='thresh',
                    default=5,
                    type=int,
                    help='Specifies the trigger threshold in units of \
                          sigma = standard deviations.')

parser.add_argument('--separ',
                    dest='separation_pulses',
                    default=100,
                    type=int,
                    help='Specifies the required separation between \
                          two pulses in a trace in units of samples.')                    


parse_args = parser.parse_args()

path_to_data_dir  = parse_args.path_to_data_dir
file_spec         = parse_args.file_spec
empty_ch          = parse_args.empty_ch
thresh            = parse_args.thresh
separation_pulses  = parse_args.separation_pulses


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
npz_dir           = '/sps/grand/pcorrea/{}/{}/{}/transient_search/'.format(array,month,mode)
npz_file_template = 'transient_search_{}_{}_{}_du_{}_thresh_{}_separ_{}.npz'.format(array,month,mode,'{}',thresh,separation_pulses)

# Get the list of used DUs
du_list = tadc.get_list_of_all_used_dus()

msg = '\nFiles contain data from following DUs: {}'.format(du_list)
print(msg)


# def get_window(trace,
#                entry,
#                width=100):

#     # Get the window of specified width around entry of interest
#     window = [entry-width/2,entry+width/2+1]

#     # If window is out of bounds of trace, set window edges to
#     # trace bounds (will reduce width)
#     bounds_trace = [0,len(trace)]

#     if window[0] < bounds_trace[0]:
#         window[0] = bounds_trace[0]

#     if window[1] > bounds_trace[1]:
#         window[1] = bounds_trace[1]

#     # The output is the trace entries corresponding to the window
#     return np.arange(window[0],window[1],dtype=int)


def search_transients(trace,
                      thresh=5,
                      separation_pulses=100):
    
    window    = []
    crossings = []

    std = np.std(trace)
    
    # Find the trace positions where N sigma is exceeded
    exceed_thresh        = np.abs(trace) > thresh*std
    samples_above_thresh = np.where(exceed_thresh)[0]

    # Stop here if there are no transients
    if not np.any(exceed_thresh):
        return window, crossings

    # Find the separation between consecutive crossings
    separation_samples = np.diff(samples_above_thresh)

    # Locate where two transients are separated
    where_separated = np.where(separation_samples > separation_pulses)[0]

    # Get the beginning of the first pulse
    begin_pulse = np.max( [0,samples_above_thresh[0]-separation_pulses/2] )
    
    # If there are multiple pulses, add different pulses
    if len(where_separated) > 0:
        for idx in where_separated:
            # Get the end of current pulse
            end_pulse = np.min( [len(trace)-1,samples_above_thresh[idx]+separation_pulses/2] )
            
            # Get window
            new_window = np.array( [begin_pulse,end_pulse], dtype=int )

            # Count how many crossings (in units of std = sigma) in window
            trace_window         = trace[new_window[0]:new_window[1]+1]
            exceed_thresh_window = np.abs(trace_window) > thresh*std
            crossings_window     = np.abs(trace_window[exceed_thresh_window]) / std

            # Add information to lists
            window.append(list(new_window))
            crossings.append(list(crossings_window))

            # Find beginning of next pulse
            begin_pulse = end_pulse - separation_pulses + separation_samples[idx]

    # Get the end of the last pulse
    end_pulse = np.min( [len(trace)-1,samples_above_thresh[-1]+separation_pulses/2] )

    # Get last window
    new_window = np.array( [begin_pulse,end_pulse], dtype=int )

    # Count how many crossings (in units of std = sigma) in last window
    trace_window         = trace[new_window[0]:new_window[1]+1]
    exceed_thresh_window = np.abs(trace_window) > thresh*std
    crossings_window     = np.abs(trace_window[exceed_thresh_window]) / std

    # Add information to lists
    window.append(list(new_window))
    crossings.append(list(crossings_window))

    return window, crossings





#########################################
# COMPUTE AND PLOT MEAN PSD FOR EACH DU #
#########################################

# Make dictionaries, that way you don't have to loop over all files
# each time you want to analyze a different DU

# Save time windows of traces (unit = sample number)
pulse_windows_x = {}
pulse_windows_y = {}
pulse_windows_z = {}

# Save the threshold crossings of transient (unit = sigma)
crossings_x = {}
crossings_y = {}
crossings_z = {}

# Save the GPS times
gps_time = {}

for du in du_list:
    pulse_windows_x[du] = np.zeros((n_entries_tot),dtype=object)
    pulse_windows_y[du] = np.zeros((n_entries_tot),dtype=object)
    pulse_windows_z[du] = np.zeros((n_entries_tot),dtype=object)

    crossings_x[du] = np.zeros((n_entries_tot),dtype=object)
    crossings_y[du] = np.zeros((n_entries_tot),dtype=object)
    crossings_z[du] = np.zeros((n_entries_tot),dtype=object)

    gps_time[du] = np.zeros((n_entries_tot),dtype=int)


# Loop over all files in run
n_files   = len(root_files)
idx_entry = 0
for i, root_file in enumerate(root_files[:1]):
    tadc  = rt.TADC(root_file)
    trawv = rt.TRawVoltage(root_file)

    n_entries_file = tadc.get_number_of_entries()

    msg = '\n{}/{}: Looping over {} events in {}'.format(i+1,n_files,n_entries_file,root_file)
    print(msg)

    # Loop over all entries in file
    for entry in range(n_entries_file):
        tadc.get_entry(entry)
        trawv.get_entry(entry)

        # Use the DU as key for the FFT dictionaries
        # NOTE: assumes that there are no coincident events
        #       i.e. only one DU (and trace etc) per event
        du = tadc.du_id[0]
    
        # Get the traces
        new_trace_ch = np.array(tadc.trace_ch[0],dtype=object)
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

        trace_x = new_trace_ch[0]
        trace_y = new_trace_ch[1]
        trace_z = new_trace_ch[2]

        new_windows_x, new_crossings_x = search_transients(trace_x,
                                                           thresh=thresh,
                                                           separation_pulses=separation_pulses)
        new_windows_y, new_crossings_y = search_transients(trace_y,
                                                           thresh=thresh,
                                                           separation_pulses=separation_pulses)
        new_windows_z, new_crossings_z = search_transients(trace_z,
                                                           thresh=thresh,
                                                           separation_pulses=separation_pulses)
        
        pulse_windows_x[du][idx_entry] = new_windows_x
        pulse_windows_y[du][idx_entry] = new_windows_y
        pulse_windows_z[du][idx_entry] = new_windows_z

        crossings_x[du][idx_entry] = new_crossings_x
        crossings_y[du][idx_entry] = new_crossings_y
        crossings_z[du][idx_entry] = new_crossings_z

        # print(entry,'Pulse windows X:',new_windows_x)
        # print(entry,'Pulse windows Y:',new_windows_y)
        # print(entry,'Pulse windows Z:',new_windows_z)

        # print(entry,'Threshold crossings X:',new_crossings_x)
        # print(entry,'Threshold crossings Y:',new_crossings_y)
        # print(entry,'Threshold crossings Z:',new_crossings_z)

        # if len(new_crossings_x) == 1:
        #     if (new_windows_x[0][1] - new_windows_x[0][0] < 120) and (new_windows_x[0][1] - new_windows_x[0][0] > 100):
        #         if np.any(np.array(new_crossings_x[0]) > 5):
        #             print(entry,new_crossings_x[0])

        # In GP13 data of May 2023, GPS of DU 1078 does not work
        if du == 1078:
            if entry == n_entries_file-1:
                trawv.get_entry(entry-1)
            else:
                trawv.get_entry(entry+1)

        gps_time[du][idx_entry] = trawv.gps_time[0]

        idx_entry += 1


for du in du_list:
    # Get rid of entries that were not filled because
    # of bad quality data, like overheating etc
    mask_filled = np.where(pulse_windows_x[du] != 0)

    pulse_windows_x[du] = pulse_windows_x[du][mask_filled]
    pulse_windows_y[du] = pulse_windows_y[du][mask_filled]
    pulse_windows_z[du] = pulse_windows_z[du][mask_filled]

    crossings_x[du] = crossings_x[du][mask_filled]
    crossings_y[du] = crossings_y[du][mask_filled]
    crossings_z[du] = crossings_z[du][mask_filled]

    gps_time[du] = gps_time[du][mask_filled]

   # Save total power in npz file
    npz_file = npz_file_template.format(du)
    np.savez(npz_dir+npz_file,
             du=du,
             pulse_windows_x=pulse_windows_x[du],
             pulse_windows_y=pulse_windows_y[du],
             pulse_windows_z=pulse_windows_z[du],
             crossings_x=crossings_x[du],
             crossings_y=crossings_y[du],
             crossings_z=crossings_z[du],
             gps_time=gps_time[du])
    msg = '\nSaved npz file: {}'.format(npz_dir+npz_file)
    print(msg)