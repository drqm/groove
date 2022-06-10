from mne_epoching import mne_epoching
import numpy as np

def multitaper(sbj, epochs):
	'''
	    This function computes a time frequency decomposition of electrophysiological, in such a way
	    to obtain high frequency broadband power(HF BB). We take a series of band ranges from 70-160. 

	    Parameters:
	        sbj (string) : Subject Identifier
	        epochs (Epoch) : MNE Epoch object

	    Returns:
	        tfr_power (Power) : Time Frequency Decomposition of epoched data, without averaging.
	                            Data in this class will be an array of following dimensions
	                            (n_epochs, n_electrodes, n_freq, n_timepoints)
	'''

	decim_parameter = 2 # Factor for decimation
	band_range = np.arange(70,160,10)
	n_cycles_inst = band_range/2
	time_bandwidth_inst = 7

	tfr_power = mne.time_frequency.tfr_multitaper(epochs, band_range, n_cycles_inst,
												  time_bandwidth=time_bandwidth_inst, 
												  use_fft=True, return_itc=False, 
												  decim=decim_parameter, average=False,
												  verbose=None, n_jobs=1)

	if not (os.getcwd() == f'/home/knight/groove/data/{sbj}'):
		os.chdir(f'/home/knight/groove/data/{sbj}')

	tfr_power.save(f"{sbj}_bothblocks_multitaper_HFA_decomposition-tfr.h5", overwrite=True)
	return tfr_power


def log_normalize(tfr_power):
	'''
	    This function takes time frequency decomposition of epoched data and first log normalizes according to 
	    a baseline of the first second of data. Then the data is averaged across frequency bands to obtain an
	    estimate of high frequency broadband (HF BB) power. NOTE: This script was written with compute limitations
	    in mind (memory). As such, it might not be as 'pythonic' as it could be.

	    Parameters:
	        tfr_power (Power) : Time Frequency Decomposition of epoched data, without averaging.
	                            Data in this class will be an array of following dimensions
	                            (n_epochs, n_electrodes, n_freq, n_timepoints)

	    Returns:
	        HFA_power_normalized (Array) : Array containing HF BB power, trial and electrode wise.
	                                       Array is of the following dimensions
	                                       (n_epochs, n_electrodes, n_timepoints)
	'''
	# tfr_power = tfr_power[0]
	(n_trials, n_electrodes, n_freq, n_timepoints) = tfr_power.data.shape
	HFA_power_normalized = np.zeros((n_trials, n_electrodes, n_timepoints))
	sfreq = int(tfr_power.info['sfreq'])
	for i in range(n_freq):
		mean_baseline = np.mean(np.log(tfr_power.data[:,:,i,0:sfreq]),axis=(0,2), keepdims=True)
		std = np.std(np.log(tfr_power.data[:,:,i,0:sfreq]), axis=(0,2), keepdims=True)
		power = (np.log(tfr_power.data[:,:,i,:])-mean_baseline) / std
		HFA_power_normalized += power

	HFA_power_normalized /= n_freq

	return HFA_power_normalized