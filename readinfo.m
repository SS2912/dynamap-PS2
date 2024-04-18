function info = readinfo(subj_code, chan_names, tabledir, varargin)

% inputs
%   subj_code: bids code of subject to match for adding info
%   chan_names: list of electrodes of that subject
%   tabledir: path of table with info on channels and subjects (for now: roi, anat, group)
% varargin
%   simEFdir: path of table with simEF values per each channel (obtained via EFtable.m function)


% % debug and test: 
% tabledir = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\ChanInfo_PS2.xlsx";
% simEFdir = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\EFtable2023-04-28_12.00.csv";
% subj_code = "2305a0fd66e8";
% load('\\dynaserv\Galvani_ps2\analysis\delphos\all_mat-thr80\baseline\sub-2305a0fd66e8_baseline.mat');
% chan_names = results.labels;

    for ii = 1:2:nargin-3
            if strcmp('simEFdir', varargin{ii})
                simEFdir = varargin{ii+1}; 
            end
    end

    chan_info = readtable(tabledir, 'Sheet', "chaninfo");
    %Associate each channel to its ROI (or other info) from the info_table
    info = [];
    subchan = chan_info(chan_info.SubjBIDS == string(subj_code), :); %extract subtable of subject's data from the excel sheet
    subchan.Channel = extractBefore(string(subchan.Channel), '-'); %only takes first electrode of the bipolar couple cause psd output is made this way
    group = string(subchan.group(1));
    % for each channel in 'channels', add info present on the excel table
    for i=1:length(chan_names)
       check = extractBefore(string(chan_names(i)), "-"); % extracts info from monopolar channel list since montages can be diff in chaninfo table and in analysed file (mono vs bipolar montage)
       
       if ~isempty(string(subchan.group(subchan.Channel == check)))
           roi = string(subchan.ROI(subchan.Channel == check));
           EImax = subchan.EImax(subchan.Channel == check);
           anat = string(subchan.Structure(subchan.Channel == check));
           
           info = [info; subj_code, chan_names(i), roi, EImax anat, group];
       else
           info = [info; subj_code, chan_names(i), "", "", "", group];
       end    
    end

    info = array2table(info, 'VariableNames',{'subj','channel','roi', 'EImax', 'anat', 'group'});
    clearvars check

    %% read simEF values and add them to other informations
    if exist('simEFdir','var') 
        simEF_table = readtable(simEFdir);
        channels_mono = extractBefore(string(info.channel), '-');
        subchan = simEF_table(simEF_table.subj == subj_code, :); 
        for i=1:length(channels_mono)
           check = string(channels_mono(i));
           simEF = string(subchan.simEF(subchan.channel == check));
           matter = string(subchan.matter(subchan.channel == check));
    
           if ~isempty(simEF)
               info.simEF(i) = simEF;
               info.matter(i) = matter;
           else
               info.simEF(i) = "";
               info.matter(i) = "";
           end    
        end

    end

end