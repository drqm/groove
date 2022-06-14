project_dir = '/home/knight/groove'

import pandas as pd
import os.path as op
scripts_path = op.join(project_dir, 'scripts')
from sys import path, argv; path.append(scripts_path)
from src.bci_processing import bci2bids

bids_root = op.join(project_dir,'data','BIDS')
SBJ = 'SLCH002'
if len(argv) > 1:
    SBJ = argv[1]

ch_info_file = op.join(project_dir,'misc/ch_info/',SBJ,'channel_types.tsv')
path_info_file = op.join(project_dir,'misc/path_info.csv')
path_info = pd.read_csv(path_info_file)

sub_ix = path_info.subject==SBJ
elec_path = path_info.elec_path.loc[sub_ix][0]
fs_subjects_dir = path_info.fs_subjects_dir.loc[sub_ix][0]
fs_subject = path_info.fs_subject.loc[sub_ix][0]

for b in [1,2]:
    bci_file = path_info['bci_file'+str(b)].loc[sub_ix][0]
    bci2bids(SBJ=SBJ, bids_root=bids_root, bci_file=bci_file,
             elec_path=elec_path, fs_subjects_dir=fs_subjects_dir,
             fs_subject=fs_subject,ch_info_file=ch_info_file,
             event_selection=list(range(1,29)),
             response_selection = [50,51,52,53,54,96,97,98,99,100],
             overwrite=True)

