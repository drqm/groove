%% 1. add Fieldtrip to the path and initialize it
ftDir = '~/fieldtrip/';
if exist('ft_defaults.m', 'file') == 0
    addpath(ftDir); ft_defaults;
end
addpath(genpath('/home/knight/groove/scripts'));

%% 2. define parameters
pos1 = [1 45 1144 1200]; % to resize ft_databrowser
pos2 = [1145 45 1144 2400];

%% 3. get subjects and specify blocks
block = 1; 
sub_code = 1;

[SBJs, ~] = get_subjects();
SBJ = SBJs{sub_code};

datadir = sprintf('/home/knight/WashU/data/BerkeleyGrooveTask/%s/BerkeleyGrooveTask/ECOG%03d/',SBJ,block);
datafiles = dir([datadir,'*.dat']);
%% 4. load data
filename = fullfile(datadir,datafiles(1).name);
hdr = ft_read_header(filename);
fsOrig = hdr.Fs;

cfg = [];
cfg.dataset = filename;
data = ft_preprocessing(cfg);
%% 5. Resample data if necessary
if fsOrig > 1000
    cfg = [];
    cfg.resamplefs = 1000;
    data = ft_resampledata(cfg, data);
end

%% 6. observe raw data
cfg = [];
cfg.viewmode = 'vertical';
cfg.preproc.demean = 'yes';
cfg.blocksize = 10;
cfg.position = pos1;
cfg.ylim = [-1 1]*50;
cfg.channel = data.label(1:50);
ft_databrowser(cfg, data);