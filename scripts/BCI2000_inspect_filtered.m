%% 1. add Fieldtrip to the path and initialize it
ftDir = '~/fieldtrip/';
if exist('ft_defaults.m', 'file') == 0
    addpath(ftDir); ft_defaults;
end
addpath(genpath('/home/knight/groove/scripts'));

%% 2. define parameters
HPcutoff = 1; % high-pass filter cutoff
LPcutoff = 250; % low-pass filter cutoff
pos1 = [1 45 1144 1200]; % to resize ft_databrowser
pos2 = [1145 45 1144 2400];

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

%% 3. load data
filename = fullfile(datadir,datafiles(1).name);
hdr = ft_read_header(filename);
fsOrig = hdr.Fs;

cfg = [];
cfg.dataset = filename;
%cfg.trl = trl;
data = ft_preprocessing(cfg);
%%
if fsOrig > 1000
    cfg = [];
    cfg.resamplefs = 1000;
    data = ft_resampledata(cfg, data);
end

%% 4. parse data
cfg = [];
cfg.channel = chantypeix{strcmp(chantypes,'iEEG')};
data = ft_selectdata(cfg, data);

%% 5. swap raw labels for anatomical labels
load(modelname);
electrodeNames = cellfun(@(x) split(x, '-'), electrodeNames, 'UniformOutput', false);
electrodeNames = cellfun(@(x) [x{1}(3:end) '_' x{end}], electrodeNames, 'UniformOutput', false);
data.label = electrodeNames;

%% 6. filter out power line noise and detrend
powerLineF0 = 60;
cfg = [];
cfg.bsfilter = 'yes';
cfg.bsfreq = (powerLineF0 * (1:round((LPcutoff+powerLineF0)/powerLineF0)) + [-1; 1])';
cfg.bsfiltord = 3;
cfg.hpfilter = 'yes';
cfg.hpfreq = HPcutoff;
cfg.hpfiltord = 3;
cfg.lpfilter = 'yes';
cfg.lpfreq = LPcutoff;
cfg.lpfiltord = 3;
data = ft_preprocessing(cfg, data);

%% 6.5. observe filtered data
cfg = [];
cfg.viewmode = 'vertical';
cfg.preproc.demean = 'yes';
cfg.blocksize = 16;
cfg.position = pos1;
cfg.ylim = [-1 1]*50;
cfg.channel = data.label(1:100);
ft_databrowser(cfg, data);

