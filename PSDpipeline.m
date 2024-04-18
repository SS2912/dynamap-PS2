function subject = PSD_PS2(filedir, tabledir, window, n_sections, varargin)

% function to calculate mean PSD and/or comparison z-value for each channel and for each subject and
% for each section (baseline, shame etc) from the raw output of the
% spectral plugin in anywave. Inputs:

% - filedir: dir of PSD .mat files with sub names at beginning and organised in sections: baseline, postA, postB, sham, postALL (tacsValidation), stimA, stimB (n_sections)
% - tabledir: directory of channel/subject information table (containing the name of the file and the extension)
% - window: window used by the spectral plugin
% - n_sections: number of sections (baseline, sham, stim, ...) 
% - varargin: baseline (default= baseline, otherwise "sham"), fmax (default=80Hz), freq_range (default= [1,45]), whichnorm (default="z-value")

for ii = 1:2:nargin-4
        if strcmp('fmax', varargin{ii})
            fmax = varargin{ii+1}; 
        elseif strcmp('freq_range', varargin{ii})
            freq_range = varargin{ii+1};
        elseif strcmp('baseline', varargin{ii})
            baseline = varargin{ii+1};
        elseif strcmp('whichnorm', varargin{ii})
            whichnorm = varargin{ii+1};
        end
end

%% debug & example:
% filedir = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\PSD\monopolar_rereferenced";
% tabledir = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\ChanInfo_PS2.xlsx";
% window = 5; n_sections = 7;

if ~exist("fmax", 'var') 
    fmax   = 80; % max freq at which to cut data for visualization and analysis
end
if ~exist("whichnorm", 'var')
    whichnorm = "z-value";
end
if ~exist("freq_range", 'var')
    freq_range = [1,45];
end
if ~exist("baseline", 'var')
    baseline = "";
end


cd(filedir) 

to_read = dir('*.mat');  % ord
chan_info = readtable(tabledir, 'Sheet', "chaninfo");
% subj_info = readtable(tabledir, 'Sheet', "Allsubj_BIDSimport");
subj_iter = 1;  %subject number for cycle (needed for the 3d matrices
subject   = struct();
psd=[];

% load psd values for each subject and calculate the difference/ratio post vs pre
for subj = 1:n_sections:length(to_read)
    subj_code = string(extractBetween(to_read(subj).name, "sub-", "_ses")); 
    
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
    
    %3a. calculate the mean PSD per each channel on the window length, for a specific frequency range
    for chan=1:length(chan_names)
        mean_baseline(chan) = mean(baseline{1, 1}(chan).psd(Frange_min:Frange_max)); % vector of fft values for that frequency f, containing all iteration (window) values
        mean_sham(chan)     = mean(sham{1, 1}(chan).psd(Frange_min:Frange_max));
        mean_stimA(chan)    = mean(stimA{1, 1}(chan).psd(Frange_min:Frange_max));
        mean_stimB(chan)    = mean(stimB{1, 1}(chan).psd(Frange_min:Frange_max));
        mean_afterA(chan)   = mean(postA{1, 1}(chan).psd(Frange_min:Frange_max));
        mean_afterB(chan)   = mean(postB{1, 1}(chan).psd(Frange_min:Frange_max));
        mean_afterAll(chan) = mean(postAll{1, 1}(chan).psd(Frange_min:Frange_max));
    end
            
            
    %3b. test before vs after via psd.fft_iterations: for each subject(iter),take channel per channel and for each frequency, test fft_iteration before vs after T
        
    for chan=1:length(chan_names)
        for f=1:samplemax

            iter_base = baseline{1, 1}(chan).fft_iterations(f,:); % vector of fft values for that frequency f, containing all iteration (window) values
            iter_sham = sham{1, 1}(chan).fft_iterations(f,:);
            iter_stimA = stimA{1, 1}(chan).fft_iterations(f,:);
            iter_stimB = stimB{1, 1}(chan).fft_iterations(f,:);
            iter_afterA = postA{1, 1}(chan).fft_iterations(f,:);
            iter_afterB = postB{1, 1}(chan).fft_iterations(f,:);
            iter_afterAll = postAll{1, 1}(chan).fft_iterations(f,:);
            
            [psham,hsham,statssham] = ranksum(iter_sham, iter_base, 'method','approximate'); % returns 1 zval for each freq, for each channel and for each patient
           
            if strcmp(baseline, "sham")
                clearvars iter_base
                iter_base = iter_sham;
            end

            [pstimA,hstimA,statsstimA] = ranksum(iter_stimA, iter_base, 'method','approximate'); 
            [pstimB,hstimB,statsstimB] = ranksum(iter_stimB, iter_base, 'method','approximate'); 
            [pA,hA,statsA] = ranksum(iter_afterA, iter_base, 'method','approximate'); 
            [pB,hB,statsB] = ranksum(iter_afterB, iter_base, 'method','approximate'); 
            [pPostAll,hPostAll,statsPostAll] = ranksum(iter_afterAll, iter_base, 'method','approximate'); 

            zvaluessham(f, chan)= statssham.zval;             % created 3d matrix with zval of the ranksum statistics for each channel and for each frequency, inside the for loop of subjects
            zvaluesStimA(f, chan)= statsstimA.zval; 
            zvaluesStimB(f, chan)= statsstimB.zval; 
            zvaluesPostA(f, chan)= statsA.zval; 
            zvaluesPostB(f, chan)= statsB.zval; 
            zvaluesPostAll(f, chan)= statsPostAll.zval; 

        end
   end

   
    %Associate each channel to its ROI (or other info) from the info_table
    info = [];
    subchan = chan_info(chan_info.SubjBIDS == string(subj_code), :); %extract subtable of subject's data from the excel sheet
    subchan.Channel = extractBefore(string(subchan.Channel), '-'); %only takes first electrode of the bipolar couple cause psd output is made this way
    % for each channel in 'channels', add info present on the excel table
    for i=1:length(chan_names)
       check = extractBefore(string(chan_names(i)), "-"); % extracts info from monopolar channel list since montages can be diff in chaninfo table and in analysed file (mono vs bipolar montage)
       
       if ~isempty(string(subchan.group(subchan.Channel == check)))
           info = [info; subj_code, check, string(subchan.ROI(subchan.Channel == check)), ...
               string(subchan.Structure(subchan.Channel == check)), string(subchan.group(subchan.Channel == check))];
       else
           info = [info; subj_code, check, "", "", ""];
       end    
    end


%    
    %create big structure with all fields
    subject(subj_iter).code     = subj_code;
    subject(subj_iter).chan     = chan_names;
    subject(subj_iter).norm     = whichnorm;
    subject(subj_iter).group    = info(1,5);
    subject(subj_iter).chaninfo = info(:,:);
    
    subject(subj_iter).zvalSHAM = zvaluessham';
    subject(subj_iter).zvalSTIMA = zvaluesStimA';
    subject(subj_iter).zvalSTIMB = zvaluesStimB';
    subject(subj_iter).zvalPOSTA = zvaluesPostA';
    subject(subj_iter).zvalPOSTB = zvaluesPostB';
    subject(subj_iter).zvalPOSTALL = zvaluesPostAll';
    
    subject(subj_iter).meanBASE = mean_baseline';
    subject(subj_iter).meanSHAM = mean_sham';
    subject(subj_iter).meanSTIMA = mean_stimA';
    subject(subj_iter).meanSTIMB = mean_stimB';
    subject(subj_iter).meanPOSTA = mean_afterA';
    subject(subj_iter).meanPOSTB = mean_afterB';
    subject(subj_iter).meanPOSTALL = mean_afterAll';
    
    subject(subj_iter).filedir = filedir;
    

    subj_iter = subj_iter+1;
    clearvars psd fft_overlap fft_window windowing subj_code f zvaluessham zvaluesStimA zvaluesStimB zvaluesPostA zvaluesPostB zvaluesPostAll ...
        mean_baseline mean_sham mean_stimA mean_stimB mean_afterA mean_afterB mean_afterAll
     
end



% %% separation in frequency bands (needs the trasnformation from bits to
% % freq, 1H/0.2=5 Hz/bits --> if data have 125 elements, it means 25 Hz (125/5) 
% avg = "median";
% 
% 
% 
% for subj = 1:length(subject)
% 
%     deltaA  = median(subject(subj).zvalSTIMA(1*window:4*window,:));
%     thetaA  = median(subject(subj).zvalSTIMA(4*window:8*window,:));
%     alphaA  = median(subject(subj).zvalSTIMA(8*window:15*window,:));
%     betaA   = median(subject(subj).zvalSTIMA(15*window:30*window,:));
%     lowgamA = median(subject(subj).zvalSTIMA(30*window:45*window,:));
%     broadA  = median(subject(subj).zvalSTIMA(1*window:45*window,:));
% 
%     deltaB  = median(subject(subj).zvalSTIMB(1*window:4*window,:));
%     thetaB  = median(subject(subj).zvalSTIMB(4*window:8*window,:));
%     alphaB  = median(subject(subj).zvalSTIMB(8*window:15*window,:));
%     betaB   = median(subject(subj).zvalSTIMB(15*window:30*window,:));
%     lowgamB = median(subject(subj).zvalSTIMB(30*window:45*window,:));
%     broadB  = median(subject(subj).zvalSTIMB(1*window:45*window,:));
% 
% 
% 
%     % now I have a matrix for each patient, for each channel with median
%     % (or mean) of zvalues within freq band
%     subject(subj).subbands_zA = [broadA; deltaA; thetaA; alphaA; betaA; lowgamA];
%     subject(subj).subbands_zB = [broadB; deltaB; thetaB; alphaB; betaB; lowgamB];
% 
%     clear delta* theta* beta* alpha* lowgam*
% end 



end