#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 07:24:51 2022

@author: dquiroga
"""
import nibabel
import mne
import os
import numpy as np

# Define paths and other variables
sub = 'SLCH002' 
surfname = 'pial'

# sub = 'BJH017' 
# surfname = 'pial.T1'

elec_path = '/home/knight/ecog/DATA_FOLDER/WashU/data/_IMAGING/' + sub + '/Electrodes'
this_subject_dir = '/home/knight/ecog/DATA_FOLDER/WashU/data/_IMAGING/' + sub
fs_sub_alias = 'segmentation'
avg_subjects_dir = '/usr/local/freesurfer_x86_64-5.3.0-modified-perl/subjects'
sig_elec_file = '/home/knight/dquiroga/groove/results/SLCH002/filtered/both_blocks/tmax16.0/HFA/multitaper/regression_linear_separated ratings/regression_timecourses/significant_electrodes.npy'

# Load the subject's MRI to obtain relevant transforms
mri_path = os.path.join(this_subject_dir,fs_sub_alias,'mri','T1.mgz')
t1 = nibabel.load(mri_path)

#get necessary transforms
ras_to_vox = mne.transforms.Transform(fro='ras', to='mri_voxel',trans=np.linalg.inv(t1.affine))
vox_to_mri = mne.transforms.Transform(fro='mri_voxel', to='mri',trans=t1.header.get_vox2ras_tkr())
mri_to_mni = mne.read_talxfm(fs_sub_alias, this_subject_dir)

# Load electrode coordinates, transform to voxel (ras) coordinates and then to MRI (ras) coordinates
elecs = os.listdir(elec_path)
ch_coords={}
for e in elecs:
    f = open(elec_path + '/' + e, "r")
    lines = f.readlines()
    lines = lines[0:-3]
    for lidx, l in enumerate(lines):
        if l != '\n':
            #transform to voxel coordinates
            ccoords = np.round(mne.transforms.apply_trans(ras_to_vox,np.array(l[0:-2].split(), dtype='float'))).astype(int)
            #transform to mri coordinates
            ccoords = mne.transforms.apply_trans(vox_to_mri,ccoords)
            # convert to meters and store
            ch_coords[e[0] + str(lidx + 1)] = ccoords/1000

# make a montage from electrode positions
montage = mne.channels.make_dig_montage(ch_coords, coord_frame='mri')
montage_mni = mne.channels.make_dig_montage(ch_coords, coord_frame='mri')

# transform to mni standard space
montage_mni.apply_trans(mri_to_mni)
montage_mni.apply_trans(mne.transforms.Transform(fro='mni_tal', to='mri', trans=np.eye(4)))

# Obtain transformed coordinates:
rcoord, mcoord = montage.get_positions()['ch_pos'], montage_mni.get_positions()['ch_pos']
ras_coord = np.array([rcoord[cch] for cch in rcoord])*1000
mni_coord = np.array([mcoord[cch] for cch in mcoord])*1000

# plot brain and electrodes
brain_kwargs = dict(cortex='low_contrast',alpha=.2,background='white')
brain = mne.viz.Brain(fs_sub_alias,subjects_dir=this_subject_dir,**brain_kwargs,surf=surfname)
brain.add_foci(coords = ras_coord, hemi='rh',color='red',scale_factor=0.25)   

brain_kwargs = dict(cortex='low_contrast',alpha=.2,background='white')
brain = mne.viz.Brain('fsaverage',subjects_dir=avg_subjects_dir,**brain_kwargs)
brain.add_foci(coords = mni_coord, hemi='rh',color='red',scale_factor=0.25)   

selecs = np.load(sig_elec_file)
selecs = selecs.astype(bool)
regs = ['intercept','harmony', 'rhythm', 'interaction', 'groove', 'pleasure']
colors = ['black', 'blue', 'red', 'purple', 'g', 'm']
for c in range(selecs.shape[1]):
    brain = mne.viz.Brain(fs_sub_alias,subjects_dir=this_subject_dir,**brain_kwargs, title=regs[c])
    brain.add_foci(coords = ras_coord[selecs[:,c]], hemi='rh',color=colors[c],scale_factor=0.25)  

# #Obtain transform to the original brain (should be identity for fsaverage)
# ctrans = mne.channels.compute_native_head_t(montage)
# trans = mne.transforms.Transform(fro='mri', to = 'head',trans=ctrans['trans'])

# ctrans_mni = mne.channels.compute_native_head_t(montage_mni)
# trans_mni = mne.transforms.Transform(fro='mri', to = 'head',trans=ctrans_mni['trans'])

# # create dummy evoked arrays for handling of electrode locations
# info = mne.create_info(ch_names = [ch for ch in ch_coords.keys()],ch_types='seeg',sfreq=1000)
# evkd = mne.EvokedArray(np.zeros((len(ch_coords.keys()),100)),info)
# evkd.set_montage(montage)

# info_mni = mne.create_info(ch_names = [ch for ch in ch_coords.keys()],ch_types='seeg',sfreq=1000)
# evkd_mni = mne.EvokedArray(np.zeros((len(ch_coords.keys()),100)),info_mni)
# evkd_mni.set_montage(montage_mni)

# # plot brain and electrodes
# brain_kwargs = dict(cortex='low_contrast',alpha=.2,background='white')
# brain = mne.viz.Brain(fs_sub_alias,subjects_dir=this_subject_dir,**brain_kwargs)
# brain.add_foci(coords = ras_coord, hemi='rh',color='red',scale_factor=0.25)   

# brain_kwargs = dict(cortex='low_contrast',alpha=.2,background='white')
# brain = mne.viz.Brain('fsaverage',subjects_dir=avg_subjects_dir,**brain_kwargs)
# brain.add_foci(coords = mni_coord, hemi='rh',color='red',scale_factor=0.25)   
    
