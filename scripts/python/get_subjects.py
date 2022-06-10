import numpy as np
import pandas as pd


def get_subjects():
  '''
      Helper function. Generates dataframe that contains all important information
      for a given patient. Consider turning this function into a script that saves a csv file.
      This would reduce having to rerun this every time we'd like to use this data.

      Parameters:

      Returns:
          SBJs (list) : List of strings where each element is a subject identifier
          info_df (DataFrame) : Dataframe containing electrode information
                                 for all subjects collected for task
  '''

  SBJs = ['SLCH002', 'BJH017']
  iEEG1 = np.arange(2, 114)
  iEEG2 = np.arange(116, 116+74)
  # Note that iEEG is the name is the original script,
  # but mni-python prefers ecog, so I renamed that channel ecog
  chantype = {}
  manual_drop = {}
  chantype['SLCH002'] = {'REF': np.arange(2),
                         'seeg': np.concatenate((iEEG1, iEEG2)).flatten(),
                         'EKG': np.array([114, 115]),
                         'eeg': np.arange(17)+190,
                         'DC': np.arange(16)+209,
                         'unknown': np.array([207, 208])}
  unknown1 = np.arange(4)
  unknown2 = np.arange(7)+224
  chantype['BJH017'] = {'REF': np.arange(2)+4,
                        'seeg': np.arange(218)+6,
                        'EKG': np.array([254,255]),
                        'eeg': np.arange(23)+231,
                        'DC': np.arange(15)+256,
                        'unknown': np.concatenate((unknown1,unknown2)).flatten()}
  # These channels are dropped because they have consistent epileptic activity or consistent slow wave delta activity
  manual_drop['SLCH002'] = ['rMTG_Amy13', 'rMTG_Amy14', 'rMTG_aHc1', 'rMTG_aHc2', 'rMTG_aHc3',
                            'rMTG_PHG2']
  manual_drop['BJH017'] = []
  infodf = pd.DataFrame.from_dict(chantype, orient='index')
  infodf['Manually Dropped Electrodes'] = manual_drop.values()
  infodf['onset epileptic'] = ''
  infodf['duration epileptic'] = ''
  infodf['block epileptic'] = ''
  
  # Important!!! Events at 161, and 212 for subject SLCH002 should be evaluated further, and may not need to be avoided.
  infodf.at['SLCH002', 'onset epileptic'] = {'block 1': np.array([31., 102., 146., 161., 212., 235., 262., 308.]),
                                            'block 2': np.array([])}
  infodf.at['SLCH002', 'duration epileptic'] = {'block 1': np.array([[1., 2., 2., 2., 1., 3., 2., 2.]]),
                                            'block 2': np.array([])}
  infodf.at['BJH017', 'onset epileptic'] = {'block 1': np.array([]),
                                            'block 2': np.array([])}
  infodf.at['BJH017', 'duration epileptic'] = {'block 1': np.array([]),
                                               'block 2': np.array([])}
  return SBJs, infodf