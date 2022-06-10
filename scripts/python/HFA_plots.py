import matplotlib.pyplot as plt
import numpy as np
from get_subjects import get_subjects
from mne_epoching import channel_selection


def HFA_single_electrode_plot(sbj, HFA_power, electrode_ind, sfreq):
	'''
	    This function serves to plot the high frequency broadband activity for a particular electrode.

	    Parameters:
				    sbj (String) : Subject identifier
				    HFA_power (ndarray) : Array of shape (n_trials, n_electrodes, n_timepoitns)
				    electrode_ind (int) : Index of the electrode whose HFA we'd like to plot
				    sfreq (int) : Sample rate of HFA data
	    Returns:
	            fig (Figure) : Figure containing the plots
	            ax (Axis Object) : ax variables containing what is in the plot.
	'''
	
	fig, ax = plt.subplots(figsize=(20,10))

	_, chantype = get_subjects()
	electrode_names, invert_idx = channel_selection(sbj, chantype)
	_, _, n_timepoints = HFA_power.shape
	tmax = int(n_timepoints/sfreq)
	pos = ax.imshow(HFA_power[:,electrode_ind,:], cmap='hot', aspect='auto', vmin=-3, vmax=3)

	ax.set_aspect('auto')

	plt.ylim([-0.5, HFA_power.shape[0]-0.5])

	plt.title(f"High Frequency Broadband Activity: {electrode_names[invert_idx[electrode_ind]]}", fontsize=30)

	fig.colorbar(pos, ax=ax)

	secs = np.arange(0, (tmax+1)*sfreq,sfreq)
	plt.xticks(sex, np.arange(-1, int(tmax+1),1), fontsize=20)
	plt.yticks(fontsize=20)

	plt.ylabel("Trial", fontsize=30)
	plt.xlabel("Time (s)", fontsize=30)

	return fig, ax