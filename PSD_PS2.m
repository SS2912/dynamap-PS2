function [PSDtable, subject] = PSD_PS2(filedir, tabledir, simEFdir, window, n_sections, varargin)

tic
% function to calculate mean PSD and/or comparison z-value for each channel and for each subject and
% for each section (baseline, shame etc) from the raw output of the
% spectral plugin in anywave. Inputs:

% - filedir: dir of PSD .mat files with sub names at beginning and organised in sections: baseline, postA, postB, sham, postALL (tacsValidation), stimA, stimB (n_sections)
% - tabledir: directory of channel/subject information table (containing the name of the file and the extension)
% - window: window used by the spectral plugin
% - n_sections: number of sections (baseline, sham, stim, ...) 
% - varargin: baseline (default= baseline, otherwise "sham"), fmax (default=80Hz), freq_range (default= [1,45]), whichnorm (default="z-value")

for ii = 1:2:nargin-5
        if strcmp('fmax', varargin{ii})
            fmax = varargin{ii+1}; 
        elseif strcmp('freq_range', varargin{ii})
            freq_range = varargin{ii+1};
        elseif strcmp('baseline', varargin{ii})
            base = varargin{ii+1};
        elseif strcmp('whichnorm', varargin{ii})
            whichnorm = varargin{ii+1};
        elseif strcmp('subbands', varargin{ii})
            subbands = varargin{ii+1};
        elseif strcmp('subbands_name', varargin{ii})
            subbands_names = varargin{ii+1};
              
        end
end

%% debug & example:
% filedir = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\PSD\monopolar_rereferenced";
% tabledir = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\ChanInfo_PS2.xlsx";
% window = 5; n_sections = 7;
% subbands = [1 4; 4 8; 8 15; 15 30; 30 45; 1 fmax]; %delta, theta, alpha, beta, lowgamma and broad
% subbands_names = ["delta", "theta", "alpha", "beta", "lowgamma", "broad"];

if ~exist("fmax", 'var') 
    fmax   = max(max(subbands)); % max freq at which to cut data for visualization and analysis
end
if ~exist("whichnorm", 'var')
    whichnorm = "z-value";
end
if ~exist("freq_range", 'var')
    freq_range = [1, fmax];
end
if ~exist("base", 'var')
    base = "";
end


cd(filedir) 

to_read = dir('*.mat');  % ord
% chan_info = readtable(tabledir, 'Sheet', "chaninfo");
% subj_info = readtable(tabledir, 'Sheet', "Allsubj_BIDSimport");
subj_iter = 1;  %subject number for cycle (needed for the 3d matrices
subject   = struct();
psd=[];
PSDtable = table();

% load psd values for each subject and calculate the difference/ratio post vs pre
for subj = 1:n_sections:length(to_read)
    subj_code = string(extractBetween(to_read(subj).name, "sub-", "_ses")); 
    clearvars baseline postA postB sham postAll stimA stimB data
    
    load('-mat', to_read(subj).name); baseline = {psd, fft_overlap, fft_window, windowing}; clearvars psd fft_overlap fft_window windowing
    load('-mat', to_read(subj+1).name); postA = {psd, fft_overlap, fft_window, windowing}; clearvars psd fft_overlap fft_window windowing
    load('-mat', to_read(subj+2).name); postB = {psd, fft_overlap, fft_window, windowing}; clearvars psd fft_overlap fft_window windowing
    load('-mat', to_read(subj+3).name); sham = {psd, fft_overlap, fft_window, windowing}; clearvars psd fft_overlap fft_window windowing
    load('-mat', to_read(subj+4).name); postAll = {psd, fft_overlap, fft_window, windowing}; clearvars psd fft_overlap fft_window windowing
    load('-mat', to_read(subj+5).name); stimA = {psd, fft_overlap, fft_window, windowing}; clearvars psd fft_overlap fft_window windowing
    load('-mat', to_read(subj+6).name); stimB = {psd, fft_overlap, fft_window, windowing}; 


    psd_values = zeros(length(psd(1).psd), size(psd,2), n_sections); % to preallocate space of matrix where to store psd values
    
    data = {baseline, postA, postB, sham, postAll, stimA, stimB};
    
       
   % 2. difference/ratio between sections
   % extract the psd values + chan names from the psd struct and then create a table with columns = channels, i = psd values
    
   for i = 1:n_sections 
    struc = data{i}{1,1}; %extracts the psd structure inside the loaded data
    psd_values(:,:,i) = [struc.psd]; % psd_values is a 3d matrix for the current subject, containing 1 column per cannel, 1 row per freq sample or each section (on z dim) (i=frequency, j=channels, z=section pre,posta, postb) 
    if i ==1 
        chan_names = string({struc.channel});
    end
    clearvars struc 
   end
    
    fs = double(fft_window)/window;
    step = (fs)/double(fft_window);
    real_freq = step:step:fs/2 - step; 
    real_freq = repmat(real_freq',1,length(chan_names));
    
    %cut the frequency to a max freq fmax (eg 100Hz to better see the changes in the low range)
    samplemax = round((fmax*length(real_freq)/(fs/2)),0);
    Frange_min = round((freq_range(1)*length(real_freq)/(fs/2)),0);
    Frange_max = round((freq_range(2)*length(real_freq)/(fs/2)),0);
    
    real_freq = real_freq(1:samplemax); 


   %% Associate each channel to its ROI (or other info) from the info_table

    info = readinfo(subj_code, chan_names, tabledir, 'simEFdir', simEFdir);

    %create big structure with all fields
    subject(subj_iter).code     = subj_code;
    subject(subj_iter).chan     = chan_names;
    subject(subj_iter).norm     = whichnorm;
    subject(subj_iter).group    = unique(rmmissing(info.group));
    subject(subj_iter).chaninfo = info(:,:);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% SAM's Plugin FOOOF test
%     % compare 2 psd and produce table
%     psd2 = baseline{1, 1}; %takes the psd structure inside of baseline
%     for comparison = 2:length(data)
%         clearvars psd1
%         psd1 = data{comparison}{1,1};
%         sam_table_fooof{1,comparison-1)= PS2_comparePSD(psd1, psd2, real_freq); 
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %3a. calculate the mean PSD per each channel on the window length, for a specific frequency range
%     for chan=1:length(chan_names)
%         mean_baseline(chan) = mean(baseline{1, 1}(chan).psd(Frange_min:Frange_max)); % vector of fft values for that frequency f, containing all iteration (window) values
%         mean_sham(chan)     = mean(sham{1, 1}(chan).psd(Frange_min:Frange_max));
%         mean_stimA(chan)    = mean(stimA{1, 1}(chan).psd(Frange_min:Frange_max));
%         mean_stimB(chan)    = mean(stimB{1, 1}(chan).psd(Frange_min:Frange_max));
%         mean_afterA(chan)   = mean(postA{1, 1}(chan).psd(Frange_min:Frange_max));
%         mean_afterB(chan)   = mean(postB{1, 1}(chan).psd(Frange_min:Frange_max));
%         mean_afterAll(chan) = mean(postAll{1, 1}(chan).psd(Frange_min:Frange_max));
%     end
                       
    %3b. test before vs after via psd.fft_iterations: for each subject(iter),take channel per channel and for each frequency, test fft_iteration before vs after T
          
    % separate iterations into bands (added 21/06/2023, to test)    
    %%%%% FREQ BY FREQ %%%%%%%
    for chan=1:length(chan_names)
        for f=1:samplemax

            iter_base = baseline{1, 1}(chan).fft_iterations(f,:); % vector of fft values for that frequency f, containing all iteration (window) values
            iter_sham = sham{1, 1}(chan).fft_iterations(f,:);
            iter_stimA = stimA{1, 1}(chan).fft_iterations(f,:);
            iter_stimB = stimB{1, 1}(chan).fft_iterations(f,:);
            iter_postA = postA{1, 1}(chan).fft_iterations(f,:);
            iter_postB = postB{1, 1}(chan).fft_iterations(f,:);
            iter_postAll = postAll{1, 1}(chan).fft_iterations(f,:);
           
            % calculate zscore against chosen baseline
            [psham,hsham,statssham] = ranksum(iter_sham, iter_base, 'method','approximate'); % returns 1 zval for each freq, for each channel and for each patient
           
            if strcmp(base, "sham")
                clearvars iter_base
                iter_base = iter_sham;
            end
            [pstimA,hstimA,statsstimA] = ranksum(iter_stimA, iter_base, 'method','approximate'); 
            [pstimB,hstimB,statsstimB] = ranksum(iter_stimB, iter_base, 'method','approximate'); 
            [pA,hA,statsA] = ranksum(iter_postA, iter_base, 'method','approximate'); 
            [pB,hB,statsB] = ranksum(iter_postB, iter_base, 'method','approximate'); 
            [pPostAll,hPostAll,statsPostAll] = ranksum(iter_postAll, iter_base, 'method','approximate'); 

            zvaluessham(f, chan)= statssham.zval;             % created 3d matrix with zval of the ranksum statistics for each channel and for each frequency, inside the for loop of subjects
            zvaluesStimA(f, chan)= statsstimA.zval; 
            zvaluesStimB(f, chan)= statsstimB.zval; 
            zvaluesPostA(f, chan)= statsA.zval; 
            zvaluesPostB(f, chan)= statsB.zval; 
            zvaluesPostAll(f, chan)= statsPostAll.zval; 
            
            psd_median_base(f, chan)= median(iter_base);  
            psd_median_sham(f, chan)= median(iter_sham);  
            psd_median_stimA(f, chan)= median(iter_stimA);  
            psd_median_stimB(f, chan)= median(iter_stimB);  
            psd_median_postA(f, chan)= median(iter_postA);   
            psd_median_postB(f, chan)= median(iter_postB);  
            psd_median_postAll(f, chan)= median(iter_postAll); 

        end
    end

    subject(subj_iter).zvalSHAM  = zvaluessham;
    subject(subj_iter).zvalSTIMA = zvaluesStimA;
    subject(subj_iter).zvalSTIMB = zvaluesStimB;
    subject(subj_iter).zvalPOSTA = zvaluesPostA;
    subject(subj_iter).zvalPOSTB = zvaluesPostB;
    subject(subj_iter).zvalPOSTALL = zvaluesPostAll;

    subject(subj_iter).psdBASE  = psd_median_base;
    subject(subj_iter).psdSHAM  = psd_median_sham;
    subject(subj_iter).psdSTIMA = psd_median_stimA;
    subject(subj_iter).psdSTIMB = psd_median_stimB;
    subject(subj_iter).psdPOSTA = psd_median_postA;
    subject(subj_iter).psdPOSTB = psd_median_postB;
    subject(subj_iter).psdPOSTALL = psd_median_postAll;

   clearvars statssham statsstimA statsstimB statsA statsB statsPostAll zvaluessham zvaluesStimA zvaluesStimB zvaluesPostA zvaluesPostB zvaluesPostAll 
   clearvars iter_base iter_sham iter_stimA iter_stimB iter_postA iter_postB iter_postAll
   clearvars psd_median_base psd_median_sham psd_median_stimA psd_median_stimB psd_median_postA psd_median_postB psd_median_postAll

   %%%%%  BY SUB_BANDS (median in each sub-band)
   for frange = 1:size(subbands,1) % not in original version, where I was just doing median on zvalue (wrong)
        f = subbands(frange,1)*window:subbands(frange,2)*window; % start is 1hz = sample num 5 (each sample is 0.2Hz)

        for chan=1:length(chan_names)    
            iter_base    = median(baseline{1, 1}(chan).fft_iterations(f,:)); % vector of fft values for that frequency f, containing all iteration (window) values
            iter_sham    = median(sham{1, 1}(chan).fft_iterations(f,:)); 
            iter_stimA   = median(stimA{1, 1}(chan).fft_iterations(f,:)); 
            iter_stimB   = median(stimB{1, 1}(chan).fft_iterations(f,:)); 
            iter_postA   = median(postA{1, 1}(chan).fft_iterations(f,:)); 
            iter_postB   = median(postB{1, 1}(chan).fft_iterations(f,:)); 
            iter_postAll = median(postAll{1, 1}(chan).fft_iterations(f,:)); 

            % calculate zscore against chosen baseline
            [psham,hsham,statssham] = ranksum(iter_sham, iter_base, 'method','approximate'); % returns 1 zval for each freq, for each channel and for each patient
           
            if strcmp(base, "sham")
                clearvars iter_base
                iter_base = iter_sham;
            end
            [pstimA,hstimA,statsstimA] = ranksum(iter_stimA, iter_base, 'method','approximate'); 
            [pstimB,hstimB,statsstimB] = ranksum(iter_stimB, iter_base, 'method','approximate'); 
            [pA,hA,statsA] = ranksum(iter_postA, iter_base, 'method','approximate'); 
            [pB,hB,statsB] = ranksum(iter_postB, iter_base, 'method','approximate'); 
            [pPostAll,hPostAll,statsPostAll] = ranksum(iter_postAll, iter_base, 'method','approximate'); 

            zvaluessham(chan, frange)= statssham.zval;  % 1 row per channel, 1 column per sub-band (in the order of "subbands")
            zvaluesStimA(chan, frange)= statsstimA.zval; 
            zvaluesStimB(chan, frange)= statsstimB.zval; 
            zvaluesPostA(chan, frange)= statsA.zval; 
            zvaluesPostB(chan, frange)= statsB.zval; 
            zvaluesPostAll(chan, frange)= statsPostAll.zval; 

            psd_median_base_bands(chan, frange)= median(iter_base);  
            psd_median_sham_bands(chan, frange)= median(iter_sham);  
            psd_median_stimA_bands(chan, frange)= median(iter_stimA);  
            psd_median_stimB_bands(chan, frange)= median(iter_stimB);  
            psd_median_postA_bands(chan, frange)= median(iter_postA);   
            psd_median_postB_bands(chan, frange)= median(iter_postB);  
            psd_median_postAll_bands(chan, frange)= median(iter_postAll);  
            
            clearvars iter_base iter_sham iter_stimA iter_stimB iter_postA iter_postB iter_postAll

        end
    end

    subject(subj_iter).subbands_zsham    = zvaluessham;
    subject(subj_iter).subbands_zstimA   = zvaluesStimA;
    subject(subj_iter).subbands_zstimB   = zvaluesStimB;
    subject(subj_iter).subbands_zpostA   = zvaluesPostA;
    subject(subj_iter).subbands_zpostB   = zvaluesPostB;
    subject(subj_iter).subbands_zpostAll = zvaluesPostAll;

    subject(subj_iter).psd_median_base_bands    = psd_median_base_bands;
    subject(subj_iter).psd_median_sham_bands    = psd_median_sham_bands;
    subject(subj_iter).psd_median_stimA_bands   = psd_median_stimA_bands;
    subject(subj_iter).psd_median_stimB_bands   = psd_median_stimB_bands;
    subject(subj_iter).psd_median_postA_bands   = psd_median_postA_bands;
    subject(subj_iter).psd_median_postB_bands   = psd_median_postB_bands;
    subject(subj_iter).psd_median_postAll_bands = psd_median_postAll_bands;
    
    
    subject(subj_iter).filedir = filedir;
   
    subj_iter = subj_iter+1;
    clearvars psd fft_overlap fft_window windowing subj_code f zvaluessham zvaluesStimA zvaluesStimB zvaluesPostA zvaluesPostB zvaluesPostAll ...
        statssham statsstimA statsstimB statsA statsB statsPostAll zvaluessham zvaluesStimA zvaluesStimB zvaluesPostA zvaluesPostB zvaluesPostAll 
    clearvars iter_base iter_sham iter_stimA iter_stimB iter_postA iter_postB iter_postAll
    clearvars psd_median_base_bands psd_median_sham_bands psd_median_stimA_bands psd_median_stimB_bands psd_median_postA_bands psd_median_postB_bands psd_median_postAll_bands
      
end


%% mixed model table for stats in R

PSDtable = [];
for subj = 1:length(subject)
    elem = size(subject(subj).chaninfo,1);
    b = table2cell(subject(subj).chaninfo);
    a = cellstr(reshape([b{:}],size(b)));
    PSDtable = [PSDtable; a, repelem("base", elem)', subject(subj).psd_median_base_bands];
    PSDtable = [PSDtable; a, repelem("sham", elem)', subject(subj).psd_median_sham_bands];
    PSDtable = [PSDtable; a, repelem("stimA", elem)', subject(subj).psd_median_stimA_bands];
    PSDtable = [PSDtable; a, repelem("stimB", elem)', subject(subj).psd_median_stimB_bands];
    PSDtable = [PSDtable; a, repelem("postA", elem)', subject(subj).psd_median_postA_bands];
    PSDtable = [PSDtable; a, repelem("postB", elem)', subject(subj).psd_median_postB_bands];
    PSDtable = [PSDtable; a, repelem("postAll", elem)', subject(subj).psd_median_postAll_bands];
end

PSDtable = array2table(PSDtable);
PSDtable.Properties.VariableNames = ["subject", "channel", "ROI", "EImax", "anat", "group", "simEF", "matter", "section", subbands_names];

PSDtable_zvalues = [];
for subj = 1:length(subject)
    elem = size(subject(subj).chaninfo,1);
    b = table2cell(subject(subj).chaninfo);
    a = cellstr(reshape([b{:}],size(b)));
    PSDtable_zvalues = [PSDtable_zvalues; a, repelem("zsham", elem)', subject(subj).subbands_zsham];
    PSDtable_zvalues = [PSDtable_zvalues; a, repelem("zstimA", elem)', subject(subj).subbands_zstimA];
    PSDtable_zvalues = [PSDtable_zvalues; a, repelem("zstimB", elem)', subject(subj).subbands_zstimB];
    PSDtable_zvalues = [PSDtable_zvalues; a, repelem("zpostA", elem)', subject(subj).subbands_zpostA];
    PSDtable_zvalues = [PSDtable_zvalues; a, repelem("zpostB", elem)', subject(subj).subbands_zpostB];
    PSDtable_zvalues = [PSDtable_zvalues; a, repelem("zpostAll", elem)', subject(subj).subbands_zpostAll];
end
PSDtable_zvalues = array2table(PSDtable_zvalues);
PSDtable_zvalues.Properties.VariableNames = ["subject", "channel", "ROI", "EImax", "anat", "group", "simEF", "matter", "section", subbands_names];


toc

end