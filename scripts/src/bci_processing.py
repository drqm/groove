#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import mne
import numpy as np
import os
import shutil
import pandas as pd
import nibabel
from os import path as op
from BCI2kReader import BCI2kReader as b2k
from mne_bids import BIDSPath, write_raw_bids

def bci_fetch_montage(elec_path, fs_subjects_dir, fs_subject):
    """This function obtains an electrode montage from provided
    channel information. Optimized for BCI2000 workflow.
    
    Parameters
    ----------
    elec_path: str
        Path to the folder where electrode coordinates are located.
        One file per shaft is read.
    fs_subjects_dir: str
        Path of the freesurfer subjects directory where
        the anatomy of current subject is located
    fs_subject: str
        Name of the subject in the freesurfer subjects directory.

    Returns
    ----------
    montage: DigMontage
            Obtained electrode montage with respective channel names and
            coordinates.
    """
    print('loading MRI')
    mri_path = op.join(fs_subjects_dir,fs_subject,'mri/T1.mgz')
    t1 = nibabel.load(mri_path)
    
    print('obtaining coordinate transforms')
    #get necessary transforms
    ras_to_vox = mne.transforms.Transform(fro='ras', to='mri_voxel',trans=np.linalg.inv(t1.affine))
    vox_to_mri = mne.transforms.Transform(fro='mri_voxel', to='mri',trans=t1.header.get_vox2ras_tkr())
    mri_to_mni = mne.read_talxfm(subject=fs_subject,subjects_dir=fs_subjects_dir)

    # Load electrode coordinates, transform to voxel (ras) coordinates and then to MRI (ras) coordinates
    
    print('transforming electrode coordinates from voxel (ras) to MRI (ras)')
    elecs = os.listdir(elec_path)
    ch_coords={}
    for e in elecs:
        print('loading elc info {}'.format(elec_path + '/' + e))
        f = open(elec_path + '/' + e, "r")
        lines = f.readlines()
        lines = lines[0:-3]
        lidx = 0
        for l in lines:
            if l != '\n':
                lidx = lidx+1
                #transform to voxel coordinates
                ccoords = np.round(mne.transforms.apply_trans(ras_to_vox,np.array(l[0:-2].split(), dtype='float'))).astype(int)
                #transform to mri coordinates
                ccoords = mne.transforms.apply_trans(vox_to_mri,ccoords)
                # convert to meters and store
                ch_coords[e[0] + str(lidx)] = ccoords/1000
     
    print('transforming to common MNI-TAL space')            
    # make a montage from electrode positions
    montage = mne.channels.make_dig_montage(ch_coords, coord_frame='mri')
    montage_mni = mne.channels.make_dig_montage(ch_coords, coord_frame='mri')

    # transform to mni standard space
    montage_mni.apply_trans(mri_to_mni)
    return montage

def get_events(sig, srate, selection=[]):
    """This function extracts events from BCI2000 states,
    typically StimulusCode or KeyDown.
    
    Parameters
    ----------
    sig: 1d numeric array
        signal where to find event onsets
    srate: numeric 
        sampling rate
    selection: array, list
        set of codes to look for and select.
    Returns
    ----------
    events: numeric array
         Sequence of event codes
    times: numeric array
         Sequence of event onsets   
    """
    x = np.squeeze(sig)
    oix = np.where(np.diff(x)!=0)[0]
    oix = np.array([o for o in oix if x[o] != 0])
    events = x[oix]
    times = oix / srate
    print('found the following unique codes:\n')
    print(np.unique(events))
    print('\nselecting events: {}'.format(selection))
    if any(selection):
        six = [e in selection for e in events]
        events = events[six]
        times = times[six]
    print('found {} matching events'.format(events.shape[0]))
    return events, times

def bci2bids(SBJ, bids_root, bci_file, fs_subjects_dir, fs_subject,
             ch_info_file, elec_path, run=None, event_selection = [],
             response_selection=[], overwrite=False):
    """This function converts a BCI2000 .dat file into BIDS format
    with an .edf extension for the task 'groove' and a given block.
    
    Parameters
    ----------
    SBJ: str
        Subject code
    bids_root: str
        Location of parent bids directory
    bci_file: str
        Path of the BCI file to be converted
    fs_subjects_dir: str
        Path of the freesurfer subjects directory where
        the anatomy of current subject is located
    fs_subject: str
        Name of the subject in the freesurfer subjects directory.
        This may or may not coincide with SBJ.
    ch_info_file: str
        Path to the subject-specific channel information .csv file.
        The file must contain a "name" column with electrode names,
        and a "type" column with electrode types.
    elec_path: str
        Path to the folder where electrode coordinates are located.
        One file per shaft is read. This path is input to the function
        "bci_fetch_montage".
    run: int
        Specify the current run of the task
    event_selection: array, list
        Set of codes to be selected from the "StimulusCode" BCI2000 state.
        If left blank, all events recorded will be selected. This parameter
        is input to the function "get_events".
    response_selection: array, list
        Set of codes to be selected from the "KeyDown" BCI2000 state.
        If left blank, all events recorded will be selected. This parameter
        is input to the function "get_events".
    overwrite: bool
        Whether to overwrite the current bids file.
    Returns
    ----------
    raw: mne.Raw object with the converted data.
    """
    #Load channel information
    print('loading channel info {}'.format(ch_info_file))
    ch_info = pd.read_csv(ch_info_file, delimiter='\t')
       
    # obtain sampling rate, events and responses
    print('loading sampling frequency, key presses and events')
    with b2k.BCI2kReader(bci_file) as test:
        sfreq = test.parameters['SamplingRate']
        key_press = test.states['KeyDown']
        triggers = test.states['StimulusCode']
        
    # create info
    info = mne.create_info(ch_names = list(ch_info.name), sfreq=sfreq, ch_types = list(ch_info.type))
    info['description'] = 'Data acquired and loaded with BCI2000'
    info['line_freq'] = 60
    montage = bci_fetch_montage(elec_path, fs_subjects_dir, fs_subject)
    print('created montage with channel names: ')
    print(montage.ch_names)
    info.set_montage(montage,on_missing='warn')
    
    ## Fix montage for EEG electrodes
    #info.set_montage(montage='standard_postfixed',on_missing='warn')
    print('Reading  data')
    with b2k.BCI2kReader(bci_file) as test:
        signal = test.signals
    
    raw = mne.io.RawArray(signal,info)
    
    # Get current block / task
    with b2k.BCI2kReader(bci_file) as test:
        stim = test.parameters['Stimuli'][0]
        
    if 'instructions1' in stim:
        crun = 'GrooveWantingToMove'
    elif 'instructions2' in stim:
        crun = 'GrooveLiking'
        
    print('current block is {}'.format(crun))
    # Adding annotations, triggers and responses
    print('adding triggers, annotations and responses')
    triggers, tonsets = get_events(sig=triggers, srate=sfreq,selection=event_selection)
    responses, ronsets = get_events(sig=key_press, srate=sfreq,selection=response_selection)
    raw.annotations.append(onset=0.01,duration=0, description=crun + '_start')
    raw.annotations.append(onset=raw.times[-1] - 0.01, duration=0, description=crun + '_end')
    raw.annotations.append(onset=tonsets,duration=0, description= ['SoundOnset_' + str(trig) for trig in triggers])
    raw.annotations.append(onset=ronsets,duration=0, description=['Response_' + str(resp) for resp in responses])
    
    montage = raw.get_montage()
    print('writing to bids')
    # make sure BIDS root is not erased
    #if os.path.exists(bids_root):
    #    shutil.rmtree(bids_root)
        
    bids_path = BIDSPath(subject=SBJ,root=bids_root,task=crun,run=run) 
    write_raw_bids(raw, bids_path, anonymize=dict(daysback=40000),
                   montage=montage, acpc_aligned=True, overwrite=overwrite,
                   allow_preload=True, format='EDF')
    print('done')
    return raw
