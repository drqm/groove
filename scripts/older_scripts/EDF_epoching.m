%% 1. add Fieldtrip to the path and initialize it
ftDir = '~/fieldtrip/';
if exist('ft_defaults.m', 'file') == 0
    addpath(ftDir); ft_defaults;
end
addpath(genpath('/home/knight/groove/scripts'));

%% 2. define parameters
flagFig = 0; % flag to create figures in the process
HPcutoff = 1; % high-pass filter cutoff
LPcutoff = 250; % low-pass filter cutoff
%pos1 = [1 45 1144 1200]; % to resize ft_databrowser
pos1 = [1 30 2288 900]; 
pos2 = [1145 45 1144 1200];

%% get subjects
block = 1; 
sub_code = 1;

[SBJs,chtypes] = get_subjects();
SBJ = SBJs{sub_code};

chantypeix = chtypes.ix{1}; % select iEEG channels
chantypes = chtypes.type;

datadir = sprintf('/home/knight/WashU/data/BerkeleyGrooveTask/%s/BerkeleyGrooveTask/ECOG%03d/',SBJ,block);
datafiles = dir([datadir,'*.dat']);
modelname = sprintf('/home/knight/WashU/data/_MODEL/%s',SBJ);
output_dir = sprintf('/home/knight/groove/data/%s', SBJ);
% Specify output folder
if ~isfolder(output_dir)
    mkdir(output_dir)
end

%% 3. load data
filename = fullfile(datadir,datafiles(1).name);
hdr = ft_read_header(filename);
pythonhdr_filename = sprintf("/home/knight/groove/data/%s/%s_hdr_%03d_block", SBJ, SBJ, block);
save(pythonhdr_filename, 'hdr', '-v7.3')
fsOrig = hdr.Fs;

cfg = [];
cfg.dataset = filename;
data = ft_preprocessing(cfg);
python_filename = sprintf("/home/knight/groove/data/%s/%s_raw_%03d_block", SBJ, SBJ, block);
save(python_filename, 'data', '-v7.3')

% % Get trigger onsets
tmin = -1;
tmax = 16;
trigger_ix = find(diff(hdr.orig.states.StimulusCode)~=0);
event_vals = double(hdr.orig.states.StimulusCode(trigger_ix));
trigger_ix = trigger_ix(ismember(event_vals,5:28));
event_vals = event_vals(ismember(event_vals,5:28));

trl = [trigger_ix + tmin*fsOrig,...
       trigger_ix + tmax*fsOrig,...
       ones(length(trigger_ix),1)*tmin*fsOrig,...
       event_vals];
   
   
 
% pad last trial
cfg = [];
cfg.dataset = filename;
cfg.trl = trl;
data = ft_preprocessing(cfg);
%%
if fsOrig > 1000
    cfg = [];
    cfg.resamplefs = 1000;
    data = ft_resampledata(cfg, data);
end

%% 4. parse data
cfg = [];
cfg.channel = chantypeix{strcmp(chantypes,'sEEG')};
data = ft_selectdata(cfg, data);
old_labels = data.label;

%% Re-reference
cfg = [];
cfg.reref = 'yes';
cfg.refchannel = 'all';
cfg.refmethod = 'bipolar';
cfg.groupchans = 'yes';
data = ft_preprocessing(cfg,data);

%% 5. swap raw labels for anatomical labels
load(modelname);
electrodeNames = cellfun(@(x) split(x, '-'), electrodeNames, 'UniformOutput', false);
electrodeNames = cellfun(@(x) [x{1}(3:end) '_' x{end}], electrodeNames, 'UniformOutput', false);
ref_labels = cellfun(@(x) split(x, '-'), data.label, 'UniformOutput', false);
ref_labels = cellfun(@(x) [x{1}], ref_labels, 'UniformOutput', false);

data.label = electrodeNames(ismember(old_labels, ref_labels));

%% 6. filter out power line noise and detrend
powerLineF0 = 60;
cfg = [];
% cfg.bsfilter = 'yes';
% cfg.bsfreq = (powerLineF0 * (1:round((LPcutoff+powerLineF0)/powerLineF0)) + [-1; 1])';
% cfg.bsfiltord = 3;
cfg.hpfilter = 'yes';
cfg.hpfreq = HPcutoff;
cfg.hpfiltord = 3;
cfg.lpfilter = 'yes';
cfg.lpfreq = LPcutoff;
cfg.lpfiltord = 3;
data = ft_preprocessing(cfg, data);

preprocessed_path = sprintf('/home/knight/groove/data/%s/%s_preprocessed_%03d_block', SBJ, SBJ, block);
save(preprocessed_path,'data')
% fiff_file = 'Subject1_preprocessed.fif';
% fieldtrip2fiff(fiff_file, data)
%% 6.5. observe data
cfg = [];
cfg.viewmode = 'vertical';
cfg.preproc.demean = 'yes';
cfg.blocksize = 16;
cfg.position = pos1;
cfg.ylim = [-1 1]*50;
cfg.channel = data.label(1:50);
ft_databrowser(cfg, data);

