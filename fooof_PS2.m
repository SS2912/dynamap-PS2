function fooof_PS2(dir_data, dir_info, varargin)

%% Analyze PSD .mat results (standard+FOOOF) and create a table with all results for all subjects
% 
% Syntax:  
%    fooof_table = fooof_PS2(input_folder, info_table, varargin)
%    
% Input format: 
%    1 folder per section, 1 file per subject with subj code, all subjects in the same folder (ex: path/sham/sub-0b12810f4f6d_sham.mat)
%
% Inputs:
%   input_folder    - input adress of folder containing all different sections' folders (if only 1 section--> must be the adress of that folder with.mat files inside)
%   info_table      - input full adress (with table name) of the clinical info table containing specific info for each channel
%
% Optional:
%   sections        - cell array of all input sections (==folders, default =1). Ex: {"baseline", "sham", "postA", "postB", "postall"}
%   fmax            - max frequency you want to calculte FOOF and psd for (default = 90)
%   f_range         - range of freqeuncy you want to calculate results for (default = [1 90])
%   baseline        - string of the sections you want as baseline (Ex: "sham", default = no baseline)
%   window          - size of window used to calculate psd (default = 5)
%
% Outputs:
%   psd_table       - table with one row per channel, one value per frequency bin (decibel and FOOOF: periodic, aperiodic)
%
% Required functions: 
%                   - readinfo_PTSD.m (author: S.S.)  
%
% Authors: Sara Simula (original: SM+EGa, first version: Sept 2023)
% last version: Oct 2023 (calculation AUC)

%% DEBUG
% dir_data = "\\dynaserv\Galvani_ps2\analysis\PSD\FOOOF";
% dir_info = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\ChanInfo_PS2.xlsx";
% sections = ["baseline", "sham", "tcsStimA", "tcsStimB", "postStimA", "postStimB", "postAll"];

%% Variables initialization
date = datestr(clock,'YYYY-mm-dd');
settings   = struct();  % fooof parameters : Use defaults thr to detect oscillations

% varargins
sections = "single"; % default is a string vector of length=1 (1 sections as in PTSD, only rest)
window   = 5;
fmin     = 1;
fmax     = 45;
% bands    = [1 4; 4 8; 8 15; 15 30; 30 45];
% band_name = ["delta", "theta", "alpha", "beta", "gamma"];

for ii = 1:2:nargin-2
        if strcmp('sections', varargin{ii})
            sections = varargin{ii+1}; 
        elseif strcmp('fmax', varargin{ii})
            fmax = varargin{ii+1}; 
        elseif strcmp('freq_range', varargin{ii})
            f_range = varargin{ii+1};
        elseif strcmp('baseline', varargin{ii})
            baseline = varargin{ii+1};
        elseif strcmp('window', varargin{ii})
            window = varargin{ii+1};
        elseif strcmp('bands', varargin{ii})
            bands = varargin{ii+1};
        elseif strcmp('band_name', varargin{ii})
            band_name = varargin{ii+1};
        end
end

%% Read files and compute FOOOF
% if exist("baseline", 'var')              % if we are comparing different sections, decide the baseline in advance
%     cd(strcat(dir_data, '\', baseline))
%     to_read_baseline = dir('*.mat'); 
% else 
%     cd(dir_data)                         % if we don't have multiple sections, just read the single one
%     to_read_baseline = dir('*.mat'); 
% end

n_sections = length(sections);

for sect = 1:n_sections
    fooof_table  = table();
    cd(strcat(dir_data, '\', sections(sect)))
    to_read = dir('*.mat');

    for subj = 1:length(to_read)
        subj_code = string(extractBetween(to_read(subj).name, "sub-", "_ses")); 
        
        data = load('-mat', to_read(subj).name); 

        % Compute real frequencies %
        chan_names = string({data.psd.channel});
        fs = double(data.fft_window)/window;
        step = (fs)/double(data.fft_window);
        real_freq = step:step:fs/2 - step;
        f_range  = [fmin 45.0001]; % i had to put the 0001 cause foof was not including 45 without it... not sure it will work in future

        % cut the frequency to a max freq fmax and a min freq fmin
        samplemax = find(floor(real_freq) == fmax, 1);
        samplemin = find(real_freq == fmin);
        real_freq_cut = real_freq(samplemin:samplemax);

        % read clinical info (1 line/channel) thanks to readinfo_PTSD.m %
        info = readinfo(subj_code, chan_names, dir_info); % there will be NaNs where same channel was not detected (less channels in info table wrt psd montage sometimes)
        varnames_info = info.Properties.VariableNames;

        % Loop on single channels %
        for chan = 1:length(chan_names)
            if sum(isinf(data.psd(chan).psd)) == 0
                chan_name = data.psd(chan).channel;
                rep = length(samplemin:samplemax); % length of frequency vector
                psd_dB = data.psd(chan).psd(samplemin:samplemax); % get psd value (in decibel) for one channel within the range of fmin-fmax
                if sum(strcmp(string(info.channel), chan_name))
                    info_chan = info(strcmp(string(info.channel), chan_name),:);
                else 
                    info_chan = array2table(NaN(1, size(info,2)));
                end
    
                % Compute fooof based on psd %
                mean_psd = mean(data.psd(chan).fft_iterations,2); % mean of psd values of all iterations
                res_fooof = fooof(real_freq', mean_psd, f_range, settings, true); % fooof computation
                fooof_oscilpart = res_fooof.fooofed_spectrum - res_fooof.ap_fit; % oscillation part of fooof
       
    
                % calculation of Area Under Curve (AUC) for each band %
                AUC = [];
    
                for band = 1:length(bands)
                    idx_freq_min = find(floor(real_freq_cut) == bands(band,1), 1);
                    idx_freq_max = find(floor(real_freq_cut) == bands(band,2), 1);
                    if band == length(bands)
                        x = real_freq_cut(idx_freq_min : idx_freq_max);
                        y = fooof_oscilpart(idx_freq_min : idx_freq_max)';
                    else
                        x = real_freq_cut(idx_freq_min : idx_freq_max -1);
                        y = fooof_oscilpart(idx_freq_min : idx_freq_max -1)';
                    end
                    curr_auc = trapz(x, y);
                    curr_band = band_name(band);
                    
                    nelem = length(x);
                    AUC = [AUC; repelem(curr_auc, nelem)', repelem(curr_band, nelem)'];
                end
                
                %% Create the table with the psd values for each freq (real_freq') for single channels and for each subj + section (bigger cycles up)
                psd_table_curr = table(repelem(subj_code, rep)', ...
                    repelem(sections(sect), rep)', ...
                    repelem(chan_names(chan), rep)', ...
                    real_freq_cut', mean_psd(samplemin:samplemax), psd_dB, ...
                    repelem(res_fooof.aperiodic_params(1), rep)', ...  % FOOOF offset
                    repelem(res_fooof.aperiodic_params(2), rep)', ...  % FOOOF exponential
                    res_fooof.ap_fit', res_fooof.fooofed_spectrum', ...
                    fooof_oscilpart', ...
                    AUC(:,1), AUC(:,2), ...              
                    'VariableNames', ...
                    {'sub', 'section', 'chan', 'freq', 'psd', 'psd_decibel', 'psd_offset', ...
                    'psd_exp', 'aperiodic_fit', 'fooof_spectrum', 'fooof_oscil_part', 'AUC', 'band_name'});
    
                info_toadd = array2table(string(table2cell(repmat(info_chan, rep, 1))), "VariableNames", varnames_info);
                fooof_table = [fooof_table; psd_table_curr info_toadd]; % save the data in the big table
                
                clearvars info_chan info_toadd psd_table_curr 
            end
        end
        clearvars  info rep chan_names chan_name psd_dB 

    end

nametable = strcat(dir_data, '\', date, '-PSD_FOOOF-', sections(sect), '_fmax', string(fmax), '.xlsx');
writetable(fooof_table, nametable)

end


end



