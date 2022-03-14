# This script serves to load BCI2000-ECOG data for use with MNE-python
# Assumptions:
# 1. Data is loaded into matlab first and saved as matfiles with magic parameters below

import numpy as np
import h5py
import scipy.io
import mne
from fooof import FOOOF
import re
import os
import matplotlib.pyplot as plt
import pandas as pd

def main():
    SBJ = 'BJH017'
    block = 1
    mne_fooof(SBJ,block)











def mne_fooof(SBJ, block):
    raw_filepath = f'/home/knight/groove/data/{SBJ}/{SBJ}_raw_%03d_block.mat' % block
    header_filepath = f'/home/knight/groove/data/{SBJ}/{SBJ}_hdr_%03d_block.mat' % block
    modelname = f'/home/knight/WashU/data/_MODEL/{SBJ}.mat'
    behavioral_data_path = f"/home/knight/groove/data/{SBJ}/{SBJ}_behavior_%03d_block.mat" % block
    Subject_raw = h5py.File(raw_filepath)
    header = h5py.File(header_filepath)
    SBJs, chantype = get_subjects()
    chantypes = chantype[SBJ]['seeg']
    os.chdir(f'data/{SBJ}')
    labels = []
    # print(len(Subject1_raw['data']['hdr']['label'][()][0]))
    for j in range(len(Subject_raw['data']['hdr']['label'][()][0])):
        #     print(Subject1_raw['data']['hdr']['label'][j])
        reference = Subject_raw['data']['hdr']['label'][0][j]
        double_reference = Subject_raw[reference][()]
        label = ''.join([chr(Subject_raw[reference][i][0]) for i in range(len(Subject_raw[reference][()]))])
        labels.append(label)
    labels = np.array(labels)
    new_labels = labels[chantypes]
    modeldata = scipy.io.loadmat(modelname)
    electrode_names = [modeldata['electrodeNames'][i][0][0].split('-') for i in range(len(modeldata['electrodeNames']))]
    electrode_names = [f'{electrode_names[i][0][2:]}_{electrode_names[i][-1][:]}' for i in range(len(electrode_names))]
    print(electrode_names)
    ref_labels = [new_labels[i].split('-')[0] for i in range(len(new_labels))]

    Fs = Subject_raw['data']['hdr']['Fs'][0][0]
    mne_info = mne.create_info(electrode_names, Fs,ch_types='seeg')
    trial_data = Subject_raw['data']['trial'][0][0]
    array_data = np.array(Subject_raw[trial_data])
    ecog_channel_data = array_data[:,chantypes]
    raw_data = mne.io.RawArray(ecog_channel_data.T, mne_info)
    behavioral_data = scipy.io.loadmat(behavioral_data_path)
    triggers = np.array(behavioral_data['behavioraldata'][0][0][4])

    # print(triggers.shape)
    # print(triggers)
    events = np.zeros((triggers.shape[0],3), dtype='int')
    events[:,0] = triggers[:,0]
    events[:,2] = triggers[:,3]
    # print(events)
    # Epoched_data = mne.Epochs(raw_data, events, tmax=16)
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
            ecog_data_bip_reref[i] = ecog_channel_data.T[i,:]
        else:
            ecog_data_bip_reref[i] = ecog_channel_data.T[i,:] - ecog_channel_data.T[ref_elec,:]

    h_freq = 250
    l_freq = 1

    filtered_ecog_data = mne.filter.filter_data(ecog_data_bip_reref, Fs, l_freq, h_freq)
    mne_info = mne.create_info(electrode_names, Fs,ch_types='ecog')
    filtered_raw_data = mne.io.RawArray(filtered_ecog_data, mne_info)
    tmax=17.
    Epoched_data = mne.Epochs(filtered_raw_data, events, tmax=17., baseline=(None, 0))
    sfreq = Epoched_data.info['sfreq']
    fooof_output_path = fr"/home/knight/groove/results/{SBJ}"
    # print(os.path.exists(fooof_output_path))
    if not os.path.exists(fooof_output_path):
        os.makedirs(fooof_output_path)
    fmax_inst = 35
    fmin_inst = 1
    psds, freqs = mne.time_frequency.psd_multitaper(Epoched_data, fmin=fmin_inst, fmax=fmax_inst, n_jobs=1,bandwidth=1, verbose=True)

    fooof_output_path = fr"/home/knight/groove/results/{SBJ}/block_%03d/fmax_{fmax_inst}" % block
    if not os.path.exists(fooof_output_path):
        os.makedirs(fooof_output_path)
    fooof_data = []
    for i, electrode in enumerate(electrode_names):
        fm = FOOOF()
        freq_range = [fmin_inst,fmax_inst]
        average_psds = np.mean(psds[:,i,:],axis=0)
        fm.report(freqs, average_psds, freq_range)
        fooof_plot = f'{electrode}_fooof.jpg'
        plt.savefig(os.path.join(fooof_output_path, fooof_plot))
        fooof_results_pdf = f'{electrode}_foof.pdf'
        fm.save_report(os.path.join(fooof_output_path, fooof_results_pdf))

        aperiodic_params = fm.get_params('aperiodic_params')
        peak_params = fm.get_params('peak_params')
        gaussian_params = fm.get_params('gaussian_params')
        fooof_error = fm.get_params('error')
        fooof_r2 = fm.get_params('r_squared')
        fooof_data.append([i,electrode,aperiodic_params, peak_params, gaussian_params, fooof_error, fooof_r2])
        plt.close()

    fooof_df = pd.DataFrame(fooof_data, columns=['index','Electrode Name','Aperiodic Parameters',
                                                 'Peak Parameters', 'Gaussian Parameters',
                                                 'FOOOF Error', r'$FOOOF R^{2}$'])

    fooof_output_file = os.path.join(fooof_output_path, f"{SBJ}_FOOOF.csv" % block)
    fooof_df.to_csv(fooof_output_file)


def get_subjects():
    SBJs = ['SLCH002', 'BJH017']
    iEEG1 = np.arange(112)+2
    iEEG2 = np.arange(74)+116
    # Note that iEEG is the name is the original script,
    # but mni-python prefers ecog, so I renamed that channel ecog
    chantype = {}
    chantype['SLCH002']= {'REF': np.arange(2),
                          'seeg': np.concatenate((iEEG1, iEEG2)).flatten(),
                          'EKG': np.array([114,115]),
                          'eeg': np.arange(17)+190,
                          'DC': np.arange(16)+209,
                          'unknown': np.array([207,208])}
    unknown1 = np.arange(4)
    unknown2 = np.arange(7)+224
    chantype['BJH017'] = {'REF': np.arange(2)+4,
                          'seeg': np.arange(218)+6,
                          'EKG': np.array([254,255]),
                          'eeg': np.arange(23)+231,
                          'DC': np.arange(15)+256,
                          'unknown': np.concatenate((unknown1,unknown2)).flatten()}
    return SBJs, chantype


if __name__ == "__main__":
    main()