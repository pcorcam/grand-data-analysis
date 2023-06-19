# system
import os
import glob
import argparse

# scipy
import numpy as np

# grandlib
import grand.io.root_trees as rt


###########################

'''

#######################################
# STORING TRACES AS NPZ AND TXT FILES #
#######################################

GRANDLIB BRANCH: `dev_io_root`

This script reads out a GRAND data file in ROOT format.
It plots the GPS parameters of all events in a data file,
where each detection unit (DU) is treated separately. 
This allows us to take a quick look at the GPS data,
particularly during on-site analyses.

'''


###########################
# DEFINE PARSER ARGUMENTS #
###########################

parser = argparse.ArgumentParser(description="Plot GPS parameters of ADC traces of GRAND events\
                                              in one file.")

parser.add_argument('--path',
                    dest='path_to_data_dir',
                    default='/sps/grand/data/gp13/mar2023/0dB/ROOTfiles/',
                    type=str,
                    help='Specifies the path of the directory containing the\
                          ROOT files to convert.')
parser.add_argument('--empty_ch',
                    dest='empty_ch',
                    default=0,
                    type=int,
                    help='Specifies which of the 4 ADC channels is not branched.')

parse_args = parser.parse_args()

path_to_data_dir = parse_args.path_to_data_dir
empty_ch         = parse_args.empty_ch

###################################
# READ DATA FILE AND OBTAIN TREES #
###################################

# Get the root files
root_files = sorted( glob.glob(path_to_data_dir+'*.root') )
#root_files = ['/sps/grand/data/gp13/GRANDfiles/GRAND.TEST-RAW-10s-ChanXYZ-20dB-14dus.20230521182903.027_dat.root']
n_files    = len(root_files) 

for i, root_file in enumerate(root_files[:]):
    msg = '\nReading root file {}/{}: {}'.format(i+1,n_files,root_file)
    print(msg)

    # Initiate TADC and TRawVoltage trees
    df    = rt.DataFile(root_file)
    tadc  = df.tadc
    trawv = df.trawvoltage

    # Get the number of entries (= events) in the file
    n_entries = trawv.get_number_of_entries()

    # Skip if there are no entries in file
    if n_entries == 0:
        msg  = '\nWARNING: No entries in file!'
        msg += 'Skipping file...'
        print(msg)
        continue

    # Get the number of samples in a trace
    tadc.get_entry(0)
    # BUG: `adc_samples_count_ch` does not always give what it should
    # n_samples_ch = tadc.adc_samples_count_ch[0]
    n_samples = len( tadc.trace_ch[0][0] )

    
    # Set the dimensions of the arrays to save
    traces_x  = np.zeros((n_entries,n_samples),dtype=float)
    traces_y  = np.zeros((n_entries,n_samples),dtype=float)
    traces_z  = np.zeros((n_entries,n_samples),dtype=float)
    event_ids = np.zeros(n_entries,dtype=int)


    # Get the traces and store them in the arrays
    for entry in range(n_entries):
        tadc.get_entry(entry)
        trawv.get_entry(entry)

        # Check which ADC channels were used
        # One of them has no antenna arm branched
        mask_ch = np.array(tadc.adc_enabled_channels_ch[0],dtype=bool)

        # Sometimes all four channels are enabled, but one is
        # not branched
        # For GP13 it is always channel 0 (at least in March...)
        if np.all(mask_ch):
            mask_ch[empty_ch] = False

        # Sometimes there seems to be a bug where there are less than
        # 3 enabled ADC channels
        if len( mask_ch[mask_ch] ) < 3:
            msg  = '\nWARNING for entry {}: adc_enabled_channels_ch = {}'.format(entry,mask_ch)
            mask_ch = np.ones(4,dtype=bool)
            mask_ch[empty_ch] = False
            msg += '\nFixing to {}\n'.format(mask_ch)
            print(msg)

        # Get the traces for the three active ADC channels
        # NOTE: store traces as ADC counts! That is how the trigger
        # will work in the end
        new_trace_ch = np.array(tadc.trace_ch[0],dtype=object)
        new_trace_ch = new_trace_ch[mask_ch]

        if new_trace_ch.shape != (3,n_samples):
            msg  = '\nWARNING for entry {}: new_trace_ch.shape = {} != {}'.format(entry,new_trace_ch.shape,(3,n_samples))
            msg += '\nSkipping event...\n'
            print(msg)
            continue

        # Above a GPS temperature of 55ºC the ADC is unreliable
        # Mostly seen in GP13 data of March 2023
        # Sometimes this yields corrupted GPS data that is <10ºC
        # Keeping data between 10-55ºC keeps all the good data
        # Might become obsolete with furder upgrades to the prototypes
        if (trawv.gps_temp[0] > 55) or (trawv.gps_temp[0] < 10): # [ºC]
            msg  = '\nWARNING for entry {}: gps_temp = {} not between 10 and 55 ºC'.format(entry,trawv.gps_temp[0])
            msg += '\nSkipping event...\n'
            print(msg)
            continue

        # Fill in the arrays
        traces_x[entry]  = new_trace_ch[0]
        traces_y[entry]  = new_trace_ch[1]
        traces_z[entry]  = new_trace_ch[2]
        event_ids[entry] = trawv.event_number


    # Get rid of entries that were not filled because
    # of bad quality data, like overheating etc
    mask_filled = np.any(traces_x,axis=1)
    traces_x  = traces_x[mask_filled]
    traces_y  = traces_y[mask_filled]
    traces_z  = traces_z[mask_filled]
    event_ids = event_ids[mask_filled]
    
    msg = '\nSaving {}/{} events...'.format(len(event_ids),n_entries)
    print(msg)

    # Define the directory to store the npz files
    # This will be done at the same dir level of the root directories
    # which are called `GRANDfiles` or `ROOTfiles`
    #npz_dir = path_to_data_dir.replace('/sps/grand/data/',
    #                                   '/sps/grand/pcorrea/trigger/data/')
    npz_dir = path_to_data_dir
    npz_dir = npz_dir.replace('GRANDfiles','npz').replace('ROOTfiles','npz')
    
    # Create the npz directory if it does not exist yet
    if not os.path.exists(npz_dir):
        msg  = 'Output directory does not exist yet\n'
        msg += 'Making directory: {}'.format(npz_dir)
        print(msg)
        os.mkdir(npz_dir)

    # Define the npz filename
    # Keep same name as root file, just change extension
    npz_file = root_file.split('/')[-1]
    npz_file = npz_file.replace('.root','.npz')

    # Save the file
    np.savez(npz_dir+npz_file,
             traces_x=traces_x,
             traces_y=traces_y,
             traces_z=traces_z,
             event_ids=event_ids)
    
    msg = '>>> npz file saved at {}'.format(npz_dir+npz_file)
    print(msg)


msg = 'Done!'
print(msg)

'''
#################
# END OF SCRIPT #
#################
'''