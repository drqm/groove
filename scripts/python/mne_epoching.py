import os
import pandas as pd
import h5py
import re
import numpy as np
import scipy.io
import mne


def mne_epoching(sbj, info_df, block, tmax_inst):
    '''
        This function loads in raw ECoG data, model data, and behavioral data and returns
        an MNE Epochs object to use for subsequent analysis. Importantly, this function assumes that raw
        BCI data has been processed in matlab first.

        Parameters:
            sbj (string) : Subject Identifier
            info_df (DataFrame) : Dataframe containing electrode information
                                 for all subjects collected for task
            block (int) : 1 or 2 for groove task
            tmax_inst (float):  Max time for a single trial

        Returns:
            epoched_data (Epochs object):  MNE Epochs object containing epoched data. Epoched data doesn't contain data
                                          of dropped electrodes.
            electrode_names (list) : List of strings. Each string is an electrode name and 
                                     its position in list corresponds to the row number of 
                                     the cleaned up data
            invert_idx (list) : Indices for all electrodes not dropped. Use this in combination
                                with electrode_names to get list of electrode_names that corresponds to epoched_data.
    '''
    # Magic file paths
    filepath = f'/home/knight/groove/data/{sbj}/{sbj}_raw_%03d_block.mat' % block
    behavioral_data_path = f"/home/knight/groove/data/{sbj}/{sbj}_behavior_%03d_block.mat" % block

    # Obtain the seeg electrodes from the metadata dataframe
    chantypes = info_df.at[sbj, 'seeg']

    
    # Load the mat file that has the raw data
    Subject1_raw = h5py.File(filepath)

    # Recover the electrode labels from the mat file
    labels = []
    for j in range(len(Subject1_raw['data']['hdr']['label'][()][0])):
        reference = Subject1_raw['data']['hdr']['label'][0][j]
        label = ''.join([chr(Subject1_raw[reference][i][0]) for i in range(len(Subject1_raw[reference][()]))])
        labels.append(label)
    labels = np.array(labels)
    new_labels = labels[chantypes]

    electrode_names, invert_idx = channel_selection(sbj, info_df)
    
    # Get important parameters, including  frequency and raw data
    Fs = Subject1_raw['data']['hdr']['Fs'][0][0]
    # This next line is important and can't be combined with the next line in current state! Due to limitations
    # of hd5py
    trial_data = Subject1_raw['data']['trial'][0][0]
    trial_data_array = np.array(Subject1_raw[trial_data])
    ecog_channel_data = trial_data_array[:, chantypes]

    bad_electrodes = np.array(electrode_names)[drop_electrodes_idx_tot]
    
    # Label some electrodes as bad, taken from the localization meeting or because labels are unknown
    mne_info = mne.create_info(electrode_names, Fs,
                               ch_types='seeg')
    mne_info['bads'] = list(bad_electrodes)

    # Make annotations to drop epochs later on
    # Onsets and durations are given in sec!!!, not samples!
    onset = np.ravel(info_df.at[sbj, 'onset epileptic'][f'block {block}'])
    duration = np.ravel(info_df.at[sbj, 'duration epileptic'][f'block {block}'])
    description = ['bad epileptic activity' for i in range(len(duration))]
    annotations = mne.Annotations(onset, duration, description)

    # Put annotations and raw data together into mne RawArray object
    raw_data = mne.io.RawArray(ecog_channel_data.T, mne_info)
    raw_data.set_annotations(annotations)

    # Obtain triggers for trials and more behavioral data
    behavioral_data = scipy.io.loadmat(behavioral_data_path)

    # Structure of behavioral data should be checked in matlab but is rhythm ,harmony, choices, rt, trigger
    triggers = np.array(behavioral_data['behavioraldata'][0][0][4])

    # Organizing data into format suitable for mne to use as event information
    events = np.zeros((triggers.shape[0], 3), dtype='int')
    events[:, 0] = triggers[:, 0]
    events[:, 2] = triggers[:, 3]

    # Bipolar re-referencing
    # We'll try bipolar referencing first.

    # Assumptions:
    # 1. Depth electrodes. Each shaft has contacts with name Letters#. Each shaft always start numbering at 1.
    ecog_data_bip_reref = np.zeros(ecog_channel_data.T.shape)
    for i, name in enumerate(electrode_names):
        #     print(electrode_names[i])
        string_num = re.findall('[0-9]+', name)
        #     print(string_num)
        #     print(electrode_names[0])
        index = re.search(string_num[-1], name).start()
        #     print(index)
        new_name = f"{name[0:index]}{int(string_num[-1])+1}"
        #     print(new_name)
        ref_elec = electrode_names.index(new_name) if new_name in electrode_names else None
        if ref_elec == None:
            ecog_data_bip_reref[i] = ecog_channel_data.T[i, :]
        else:
            ecog_data_bip_reref[i] = ecog_channel_data.T[i, :] - ecog_channel_data.T[ref_elec, :]

    # Band pass filtering and notch filtering to remove HF noise and power line noise
    h_freq = 250
    l_freq = 1

    filtered_ecog_data = mne.filter.filter_data(ecog_data_bip_reref, Fs, l_freq, h_freq)
    # Notch filter to remove power line noise at 60,120,180,240...
    filtered_ecog_data = mne.filter.notch_filter(filtered_ecog_data, Fs, np.arange(60, 241, 60))

    # With data re-referenced and filtered, we package everything together into an mne RawArray object
    mne_info = mne.create_info(electrode_names, Fs, ch_types='seeg')
    filtered_raw_data = mne.io.RawArray(filtered_ecog_data, mne_info)
    filtered_raw_data.set_annotations(annotations)

    channel_picks = np.array(electrode_names)[invert_idx]

    # We obtain our epoched data, using filtered_raw_data, event information, and electrode selection
    epoched_data = mne.Epochs(filtered_raw_data, events, tmax=tmax_inst,
                              baseline=(None, 0), picks=channel_picks)

    # Close h5py file (good code practice)
    Subject1_raw.close()
    return epoched_data, electrode_names, invert_idx


def merge_blocks(sbj, info_df):
    '''
        This function serves to perform epoching over one subject, but over both blocks, and then concatenate them into
        a greater epochs object.

        Parameters:
            sbj (string) : Subject Identifier
            info_df (DataFrame) : Dataframe containing electrode information
                                  for all subjects participating in task

        Returns:
            merged_epochs (Epochs) : MNE Epochs object containing electrode data
                                     from both blocks
            electrode_names (list) : List of strings. Each string is an electrode name and 
                                     its position in list corresponds to the row number of 
                                     the cleaned up data
            invert_idx (list) : Indices for all electrodes not dropped. Use this in combination
                                with electrode_names to get list of electrode_names that corresponds to epoched_data.

    '''

    if not (os.getcwd() == f'/home/knight/groove/data/{sbj}'):
        os.chdir(f'/home/knight/groove/data/{sbj}')

    block = 1
    tmax = 17.
    # Epoching
    epoched_data, electrode_names, invert_idx = mne_epoching(sbj, info_df, block, tmax)
    epoched_data_block_2, _, _ = mne_epoching(sbj, info_df, 2, tmax)

    merged_epochs = mne.concatenate_epochs([epoched_data, epoched_data_block_2], add_offset=True)
    return merged_epochs, electrode_names, invert_idx


def channel_selection(sbj, info_df):
    '''
        This functions serves to load in relevant information about subject and their electrodes and return a list of all electrode
        names and which ones are to be used in analysis

        Parameters:
            sbj (string) : Subject identifier
            info_df (DataFrame) : DataFrame containing electrode information 
                                  for all subjects participating in task

        Returns:
            electrode_names (list) : List of strings. Each string is an electrode name and 
                                     its position in list corresponds to the row number of 
                                     the cleaned up data
            invert_idx (list) : Indices for all electrodes not dropped. Use this in combination 
                                with electrode_names to get list of electrode_names that corresponds to epoched_data.
    '''

    # Magic file name
    modelname = f'/home/knight/WashU/data/_MODEL/{sbj}.mat'

    # Load in the model file and replace the electrode labels with the actual
    # electrode names
    modeldata = scipy.io.loadmat(modelname)
    electrode_names = [modeldata['electrodeNames'][i][0][0].split('-') for i in range(len(modeldata['electrodeNames']))]
    electrode_names = [f'{electrode_names[i][0][2:]}_{electrode_names[i][-1][:]}' for i in range(len(electrode_names))]
    
    # Obtain the manually dropped electrodes from the metadata dataframe
    manual_drop = np.asarray(info_df.at[sbj, 'Manually Dropped Electrodes'], dtype='object')
    
    reg_ex_unknown = re.compile('(Unknown)$')

    drop_electrodes_idx = [idx for idx,name in enumerate(modeldata['SecondaryLabel'])
                           if len(name[0]) == 1 and re.match(reg_ex_unknown,name[0][0][0][0])]
    print(f'{len(drop_electrodes_idx)} dropped electrodes due to secondary label of Unknown')
    print('These are the electrodes:')
    print(np.array(electrode_names)[drop_electrodes_idx])

    man_drop_idx = [electrode_names.index(name) for name in manual_drop]
    print(f'{len(man_drop_idx)} dropped electrodes due to epileptic or slow wave activity')
    print('These are the electrodes:')
    print(np.array(electrode_names)[man_drop_idx])

    drop_electrodes_idx_tot = list(set(man_drop_idx+drop_electrodes_idx))
    print(f'Total dropped electrodes: {len(drop_electrodes_idx_tot)}')
    
    # Channel selection
    invert_idx = [idx for idx, _ in enumerate(electrode_names) if idx not in drop_electrodes_idx_tot]
    return electrode_names, invert_idx