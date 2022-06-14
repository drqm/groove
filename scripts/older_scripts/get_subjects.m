function [SBJs, chantype] = get_subjects() 
SBJs = {'SLCH002', 'BJH017'};
chantype.type = {'REF','sEEG','EKG','EEG','DC','unknown'};
chantype.ix = {{1:2, [3:114,117:190], 115:116, 191:207, 210:225, 208:209},
               {5:6, 7:224, 255:256, 234:254, 257:272, [1:4, 225:231]}};