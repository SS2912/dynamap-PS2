%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%     PS2 pipeline      %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sara simula, Feb 2023
% pipeline of analysis for PS2 data (EF values, PSD, FC, IEDs)

%% initialization and settings
clear
clc
date = datestr(clock,'YYYY-mm-dd_HH.MM');

%% EF validation: 
% NE provides us with tables of simulated EF and we can compare them with
% the PSD calculated in the SEEG contacts at the frequency of stimulation
% of tACS (10HZ). 

% creation of channel info table (one line per channel, combination with EF values)
% dir_data = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani";
dir_data = "\\dynaserv\Galvani_ps2\simulatedEF";
EFtable = EFtable(dir_data, 'monopolar', 1);

name = strcat(['C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\' ...
    'MyProjects\Galvani\analysis\'], 'EFtable', date, '.csv');
writetable(EFtable, name)  

simEFdir = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\" + ...
    "MyProjects\Galvani\analysis\EFtable2023-10-18_13.38.csv";
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSD 
% calculation of mean PSD from spectral plugin results 
clear
simEFdir = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\" + ...
    "MyProjects\Galvani\analysis\EFtable2023-10-18_13.38.csv";
clc
date = datestr(clock,'YYYY-mm-dd_HH.MM');
filedir  = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\" + ...
    "MyProjects\Galvani\analysis\PSD\monopolar_rereferenced"; % input .mat files
tabledir = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\" + ...
    "MyProjects\Galvani\ChanInfo_PS2.xlsx";

window   = 5; 
n_sect   = 7; %now PSDpipeline func is only made for 7 sections and not less (baseline, postA, postB, sham, postAll, stimA, stimB). launch it manually if less.
baseline = "baseline";

subbands = [1 4; 4 8; 8 15; 15 30; 30 45; 1 45; 9 11; 11 15]; 
subbands_name = ["delta", "theta", "alpha", "beta", "lowgamma", "broad", "stimfreq", "highalpha"];
% fmax = 45;
fmax   = max(max(subbands));
% freq_range = [1, 45];

[PSDzvalue, PSDall]  = PSD_PS2(filedir, tabledir, simEFdir, window, n_sect, ...
    'baseline', baseline, 'subbands', subbands, 'subbands_name', subbands_name); % varargin: 'baseline' ('baseline' (default) or "sham"),  'fmax' (default=80 Hz), 'freq_range' (default= [1,45]), 'whichnorm' (default="z-value")
nametable = strcat(['C:\Users\simula\OneDrive - Aix-Marseille Université\PhD' ...
    '\MyProjects\Galvani\analysis\PSD\'], date, '_PSD_subbandsVS', baseline, '_fmax', string(fmax), '.xlsx');
namemat   = strcat(['C:\Users\simula\OneDrive - Aix-Marseille Université\PhD' ...
    '\MyProjects\Galvani\analysis\PSD\'], date, '_PSDtableVS', baseline, '_fmax', string(fmax), '.mat');
writetable(PSDzvalue, nametable)
save(namemat, "PSDall")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FOOOF
clear; clc
dir_data = "\\dynaserv\Galvani_ps2\analysis\PSD\FOOOF";
dir_info = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\ChanInfo_PS2.xlsx";
sections = ["baseline", "sham", "tcsStimA", "tcsStimB", "postStimA", "postStimB", "postAll"];
% sections = "postStimB";
bands    = [1 4; 4 8; 8 13; 13 30; 30 45]; % ONLY WORKS WITH 45 max now, otherwise change it in hard code
band_name = ["delta", "theta", "alpha", "beta", "gamma"];

fooof_PS2(dir_data, dir_info, 'sections', sections, 'window', 5, 'bands', bands, 'band_name', band_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FC 
clear
clc
close all
dir_info = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\ChanInfo_PS2.xlsx";
simEFdir = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\" + ...
    "MyProjects\Galvani\analysis\EFtable2023-10-18_13.38.csv";

thr_h2 = 0;  % thr_h2=0 for strength, >0 for degrees. Default: strength
thr_lag = 0; % value in ms, optional input. Default: 0 
norm =1; %% default, normalise strength by number of channles of each subject
roi = "all"; %"all" (default, h2 between all channels) or "EZPZ" or "NI" to calculate node strength only in subset of EZ channels (or non-inv channels)
date = datestr(clock,'YYYY-mm-dd');
% bands = ["broad", "delta", "theta", "alpha", "beta", "lowgamma"];

% calculate node strength and zvalue of node strength for sham, postA and
% postB (still need to add during stim and postAll)

% home:
% bands = "broad";
% dir_info = "C:\Users\saras\Desktop\PhD\PS2_article\ChanInfo_PS2.xlsx";
% simEFdir = "C:\Users\saras\Desktop\PhD\PS2_article\EFtable2023-10-18_13.38.csv";
% dir_data = "C:\Users\saras\Desktop\PhD\PS2_article\FC\h2\raw\broad";
% baseline = 'baseline';
% FCtable = FC_PS2(dir_data, dir_info, simEFdir, 'graph', 1, 'base', baseline, 'sorted', 1, 'instim', 0, 'roi', roi);

% broadband
% only base, sham and post stim (for tACS too):
graphs = 0; 

for band = 1:length(bands)
    band_name = bands(band);
    dir_data = strcat("C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\FC\h2\raw\", band_name);

    baseline = 'baseline';
    FCtable = FC_PS2(dir_data, dir_info, simEFdir, 'graph', graphs, 'base', baseline, ...
        'sorted', 1, 'instim', 0, 'roi', roi);
    FCtable_ezpz = FC_PS2(dir_data, dir_info, simEFdir, 'graph', graphs, 'base', baseline, ...
        'sorted', 1, 'instim', 0, 'roi', "EZPZ"); 
     clearvars name
    name = strcat('C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\FC\', 'FCtable_', band_name, '_vs_', baseline, "-", roi, "-", date, '.xlsx'); % when saving to csv, it doesnt save variable names
    name_ezpz = strcat('C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\FC\', 'FCtable_', band_name, '_vs_', baseline, "- EZPZ -", date, '.xlsx'); % when saving to csv, it doesnt save variable names
    writetable(FCtable, name)  
    writetable(FCtable_ezpz, name_ezpz) 
    
    clearvars baseline
    baseline = 'sham';
    FCtable_sham = FC_PS2(dir_data, dir_info, simEFdir, 'graph', graphs, 'base', baseline, ...
        'sorted', 1, 'instim', 0, 'roi', roi);
    
    name_sham= strcat('C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\FC\', 'FCtable_', band_name, '_vs_', baseline, '-', roi, '-', date, '.xlsx'); % when saving to csv, it doesnt save variable names
    writetable(FCtable_sham, name_sham)  

%     %also during stim
%     clearvars dir_data name
%     dir_data = strcat("C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\FC\h2\raw\WithStim\", band_name);
%     baseline = 'baseline';
%     FCtable = FC_PS2(dir_data, dir_info, simEFdir, 'graph', graphs, 'base', baseline, ...
%         'sorted', 1, 'instim', 1, 'roi', roi);
%     name = strcat('C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\FC\', 'FCtable_WITH_STIM_', band_name, '_vs_', baseline, "-", roi, "-", date, '.xlsx'); % when saving to csv, it doesnt save variable names
%     writetable(FCtable, name)  

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% graph measures DRAFT
clear
cfg = struct();% or cfg.deg_thr = 0.2; % or other number
symmetric = 1; % symmetrize (1) or not (0) the connectivity matrix to make it Undirected (WU,BU)
graphs = 0;

% baseline
h2 = load('C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\FC\h2\raw\broad\allsubj_broadband\sub-0b12810f4f6d_ses-01_task-baseline_run-01_ieeg.vhdr_algo-H2_hp-1Hz_lp-45Hz.mat');
% conn_mat = mean(h2.aw_h2,3);
conn_mat = h2.aw_h2;
baseline_top = apply_topology_measures_on_matrix(conn_mat, cfg, symmetric, graphs);

% postB
h2 = load('C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\FC\h2\raw\broad\allsubj_broadband\sub-0b12810f4f6d_ses-01_task-postStimB_run-01_ieeg.vhdr_algo-H2_hp-1Hz_lp-45Hz.mat');
conn_mat = h2.aw_h2;
postB_top = apply_topology_measures_on_matrix(conn_mat, cfg, symmetric, graphs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SPIKE RATE
clearvars -except simEFdir
clc
date = datestr(clock,'YYYY-mm-dd');

input_folder  = "\\dynaserv\Galvani_ps2\analysis\delphos\all_mat-thr80";
mtgs_folder = "\\dynaserv\Galvani_ps2\montages\delphos\*.mtg";
tabledir =  "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\ChanInfo_PS2.xlsx";
sections = {"baseline", "sham", "postA", "postB", "postall"};

[results_all, delphos_melttable] = delphos_PS2(input_folder, mtgs_folder, sections, 'tabledir', tabledir, 'simEFdir',simEFdir);
% delphos_melttable.Properties.VariableNames = ["subj", "chanRaw", "chanDelphos", "section", "alpha", "beta", "fastripple", "gamma", "ripple", "spike"];

name_delphos = strcat('C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\Delphos\', 'DelphosTable', date, '.xlsx');
writetable(delphos_melttable, name_delphos)  

