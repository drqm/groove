import scipy.behavioral_data
import numpy as np

def get_behavior_data(sbj, block, mne_object):
'''
    This function loads in behavioral data(mat file) and gets back trial onset times, and experimental conditions.
    This function is specific to the groove task.
    This function obtains rhythm, harmony and choice for each trial.

    Parameters:
        sbj (string) : Subject Identifier
        block (int) : Experimental Block
        mne_object (MNE Object) : Any mne object with a drop_log attribute.
'''

behavioral_data_path = f"/home/knight/groove/data/{sbj}/{sbj}_behavior_%03d_block.mat" % block
behavioral_data = scipy.io.loadmat(behavioral_data_path)
triggers = np.array(behavioral_data['behavioraldata'][0][0][4])
rhythms = np.unique(behavioral_data['behavioraldata'][0,0]['rhythm'])
harmonies = np.unique(behavioral_data['behavioraldata'][0,0]['harmony'])
choices = np.unique(behavioral_data['behavioraldata'][0,0]['choices'])
# epoched_data.drop_log.dtype

filtered_rhythms = [rhythm[0][0] for i, rhythm in 
                    enumerate(behavioral_data['behavioraldata'][0,0]['rhythm']) if not mne_object.drop_log[i]]
filtered_harmonies = [harmony[0][0] for i, harmony in 
                    enumerate(behavioral_data['behavioraldata'][0,0]['harmony']) if not mne_object.drop_log[i]]
filtered_choices = [choice[0] for i, choice in 
                    enumerate(behavioral_data['behavioraldata'][0,0]['choices']) if not mne_object.drop_log[i]]

rhythms_sorted_indices = np.argsort(filtered_rhythms)
harmonies_sorted_indices = np.argsort(filtered_harmonies)
choices_sorted_indices = np.argsort(filtered_choices)

return (rhythms_sorted_indices, harmonies_sorted_indices, choices_sorted_indices), (filtered_rhythms, filtered_choices, filtered_harmonies)