function [SBJs, chantype] = get_subjects() 
SBJs = {'SLCH002'};
chantype.type = {'REF','iEEG','EEG','EKG','DC','unknown'};
chantype.ix = {{1:2, [3:114,117:190], 115:116, 191:207, 210:225, 208:209}};