# system
import os
import glob
import logging
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
It computes the mean power spectrum density (PSD) of the traces
in the data file, where each detection unit (DU) is treated
separately. In particular, the PSD are correctly normalized, such
that they can be compared to other data sets/runs. The mean PSDs
are saved in npz files for each DU.

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


parse_args = parser.parse_args()

path_to_data_dir  = parse_args.path_to_data_dir


################################
# SET SOME EXPECTED PARAMETERS #
################################

n_samples = 1024
adc_range = 2**13 # [LSB]


logging.basicConfig(filename=path_to_data_dir+'data_quality.log',
                    encoding='utf-8',
                    ßlevel=logging.DEBUG)

root_files = sorted( glob.glob(path_to_data_dir+'*.root') )
n_files    = len(root_files)

for i, root_file in enumerate(root_files[10:20]):
    tadc  = rt.TADC(root_file)
    trawv = rt.TRawVoltage(root_file)
    
    n_entries_file = tadc.get_number_of_entries()

    msg = '\n{}/{}: Looping over {} events in {}'.format(i+1,n_files,n_entries_file,root_file)
    print(msg)

    warning = False

    for entry in range(n_entries_file):
        tadc.get_entry(entry)
        trawv.get_entry(entry)
        
        new_trace_ch = np.array(tadc.trace_ch[0],dtype=object)

        # Check for length of traces
        if len(new_trace_ch[0]) != n_samples:
            msg  = '\nEntry {}, DU = {}'.format(entry,tadc.du_id[0])
            msg += '\n{} samples in trace ch0 != {}'.format(len(new_trace_ch[0]),n_samples)
            msg += '\nLogging...'
            print(msg)
            logging.warning(msg)
            warning = True

        if len(new_trace_ch[1]) != n_samples:
            msg  = '\nEntry {}, DU = {}'.format(entry,tadc.du_id[0])
            msg += '\n{} samples in trace ch1 != {}'.format(len(new_trace_ch[1]),n_samples)
            msg += '\nLogging...'
            print(msg)
            logging.warning(msg)
            warning = True

        if len(new_trace_ch[2]) != n_samples:
            msg  = '\nEntry {}, DU = {}'.format(entry,tadc.du_id[0])
            msg += '\n{} samples in trace ch2 != {}'.format(len(new_trace_ch[2]),n_samples)
            msg += '\nLogging...'
            print(msg)
            logging.warning(msg)
            warning = True
        
        if len(new_trace_ch[3]) != n_samples:
            msg  = '\nEntry {}, DU = {}'.format(entry,tadc.du_id[0])
            msg += '\n{} samples in trace ch3 != {}'.format(len(new_trace_ch[3]),n_samples)
            msg += '\nLogging...'
            print(msg)
            logging.warning(msg)
            warning = True


        # Check for wrongly coded ADC values
        if np.any( np.where( np.abs(new_trace_ch[0]) > adc_range ) ):
            msg  = '\nEntry {}, DU = {}'.format(entry,tadc.du_id[0])
            msg += '\nTrace ch0 contains entries with ADC counts > {}'.format(adc_range)
            msg += '\nChecking GPS temperature: {:.3f} ºC'.format(trawv.gps_temp[0])
            msg += '\nLogging...'
            print(msg)
            logging.warning(msg)
            warning = True

        # Check for wrongly coded ADC values
        if np.any( np.where( np.abs(new_trace_ch[0]) > adc_range ) ):
            msg  = '\nEntry {}, DU = {}'.format(entry,tadc.du_id[0])
            msg += '\nTrace ch1 contains entries with ADC counts > {}'.format(adc_range)
            msg += '\nChecking GPS temperature: {:.3f} ºC'.format(trawv.gps_temp[0])
            msg += '\nLogging...'
            print(msg)
            logging.warning(msg)
            warning = True

        # Check for wrongly coded ADC values
        if np.any( np.where( np.abs(new_trace_ch[0]) > adc_range ) ):
            msg  = '\nEntry {}, DU = {}'.format(entry,tadc.du_id[0])
            msg += '\nTrace ch2 contains entries with ADC counts > {}'.format(adc_range)
            msg += '\nChecking GPS temperature: {:.3f} ºC'.format(trawv.gps_temp[0])
            msg += '\nLogging...'
            print(msg)
            logging.warning(msg)
            warning = True

        # Check for wrongly coded ADC values
        if np.any( np.where( np.abs(new_trace_ch[0]) > adc_range ) ):
            msg  = '\nEntry {}, DU = {}'.format(entry,tadc.du_id[0])
            msg += '\nTrace ch3 contains entries with ADC counts > {}'.format(adc_range)
            msg += '\nChecking GPS temperature: {:.3f} ºC'.format(trawv.gps_temp[0])
            msg += '\nLogging...'
            print(msg)
            logging.warning(msg)
            warning = True

    
    if not warning:
        msg = '\nAll data is good in file!'
        print(msg)
        logging.info(msg)


msg = '\nDONE!'
print(msg)
