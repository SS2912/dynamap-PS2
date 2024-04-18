function FCtable = FC_PS2(dir_data, dir_info, simEFdir, varargin)

%% Calculates node strength for each channel and subj in dir_data
% Syntax:  
%    FCtable = FC_PS2(dir_data, dir_info, varargin)
%
% Inputs:
%   dir_data    - input adress of folder containing h2 results from bids pipeline with sub-code at beginning and task marked (check for postA and B) 1 folder per section, 1 file per subject with subj code, all subjects in the same folder 
%   dir_info    - dir of info table with info on each channel and randomisation group (ChanInfo_PS2.xlsx)
%   simEFdir    - directory of table (with file name and extension) containing simulated EF values for each channel (obtained via EFtable)
% 
% Varargin:
%   thr_h2, thr_lag - optional thershold for h2 (thr_h2=0 for strength, >0 for degrees) and lag. Default: thr_h2 = 0; thr_lag = 0
%   base            - baseline to use for zvalue (default: base = 'baseline')
%   norm            - 1 if you want to normalise each subj by number of channels (Default= 1)
%   sorted          - outputs FC results sorting channels by simulated EF (default = 0)
%   graph           - 1 (default) to show connectivity matrices of zvalues compared to baseline chosen in base, 0 to not show any graph
%   instim          - 0 (default) to calculate only baseline, sham, postA and postB. 1 to calculate also during stim (stimA and stimB) -- only for tDCS in boradband or excluding 10Hz band in tACS subjects
%   roi             - "all" (default, h2 between all channels) or "EZ" or "NI" to calculate node strength only in subset of EZ channels (or non-inv channels)
%
% Output:
%   FCtable     - 1 row per channel (and per subject) with simEF (inputed via ChanInfo_PS2 table), mean node stength per section OUT and TOT, zvalue of TOT stenght per section respective to base 
%
% Required functions: 
%   - ins_countlinks.m
%   - ins_findmaxh2.m (contained inside countlinks)
% Authors: Sara Simula (original: Jan 2023. Last version: April 23)

% Example/debug:
% clear
% dir_data = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\FC\h2\raw\allsubj_broadband";
% dir_info = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\ChanInfo_PS2.xlsx";
% simEFdir = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\EFtable2023-04-28_12.00.csv";
% roi = "EZPZ";

% clear
% simEFdir = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\EFtable2023-04-28_12.00.csv";
% dir_info =  "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\ChanInfo_PS2_ROItest.xlsx";
% dir_data = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\FC\h2\raw\delta";
% 
% thr_h2 = 0;  % thr_h2=0 for strength, >0 for degrees. Default: strength
% thr_lag = 0; % value in ms, optional input. Default: 0 
% norm =1; %% default, normalise strength by number of channles of each subject
% base = 'sham';
% instim = 1;
% sorted = 1;
% graph = 1;
% roi = "all";


% WITH STIMA AND STIMB: dir_data = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\FC\h2\raw\broadband_onlytDCS";

% 1. Optional variables: default values
thr_h2 = 0;  % thr_h2=0 for strength, >0 for degrees. Default: strength
thr_lag = 0; % value in ms, optional input. Default: 0 % TO CHECK???
base = 'baseline';
norm =1; %% normalise strength by number of channels
sorted = 0;
graph = 0;
instim = 0;
roi = "all";

for ii = 1:2:nargin-3
        if strcmp('thr_h2', varargin{ii})
            thr_h2 = varargin{ii+1}; 
        elseif strcmp('thr_lag', varargin{ii})
            thr_lag = varargin{ii+1};
        elseif strcmp('base', varargin{ii})
            base = varargin{ii+1};
        elseif strcmp('norm', varargin{ii})
            norm = varargin{ii+1};
        elseif strcmp('sorted', varargin{ii})
            sorted = varargin{ii+1};
        elseif strcmp('graph', varargin{ii})
            graph = varargin{ii+1};
        elseif strcmp('instim', varargin{ii})
            instim = varargin{ii+1};
        elseif strcmp('roi', varargin{ii})
            roi = varargin{ii+1};
        end
end


% 2. set directory for h2 files to analyse
cd(dir_data)
myfiles = dir('*.mat');
chan_info = readtable(simEFdir);


% 3. Calculate the strength/degrees for each channel and each patient and each channel
subj_info= [];
widedata=[];

if instim
    nsec = 7;
else 
    nsec = 5;
end

for index=1:nsec:length(myfiles)
    baseline = load(myfiles(index).name); 
    postA = load(myfiles(index+1).name);
    postB = load(myfiles(index+2).name);
    sham = load(myfiles(index+3).name);
    postAll = load(myfiles(index+4).name);
    if instim
        stimA = load(myfiles(index+5).name);
        stimB = load(myfiles(index+6).name);
    end

    subj = string(extractBetween(myfiles(index).name, 'sub-', '_ses'));
    channels = string(postB.electrode_names);
   
    % read and add info on subject and channels + Simulated EF
    info = readinfo(subj, channels, dir_info, 'simEFdir', simEFdir);
    varnames = info.Properties.VariableNames;

    b = table2cell(info);
    subj_info = string(cellstr(reshape([b{:}],size(b))));

    % 3.2. Reduce matrix of h2 and lag in case "roi" is different from "all"
    if ~strcmp("all", roi)
    switch roi
        case "EZ"
            good_labels = "EZ";
        case "EZPZ"
            good_labels = '[EP.]Z'; %any character in [] followed by Z (matches "EZ" or "PZ")
        case "NI"
            good_labels = "NI";
        otherwise
            disp("incorrect roi input, possible values are: EZ, EZPZ or NI")
    end

        if instim
            H = {baseline, sham, postA, postB, postAll, stimA, stimB};
        else 
            H = {baseline, sham, postA, postB, postAll};
        end
        for dd= 1:length(H)
            idx = ~cellfun(@isempty,regexp(subj_info(:,3), good_labels));
            H{dd}.aw_h2 = H{dd}.aw_h2(idx,idx,:);
            H{dd}.aw_lag = H{dd}.aw_lag(idx,idx,:);
            H{dd}.electrode_names = H{dd}.electrode_names(idx);
        end
        if instim 
           baseline = H{1}; sham = H{2}; postA = H{3}; postB = H{4}; postAll = H{5}; stimA = H{6}; stimA = H{7};
        else
            baseline = H{1}; sham = H{2}; postA = H{3}; postB = H{4}; postAll = H{5};
        end
%         channels = string(postB.electrode_names); %update the selected channels
        subj_info = subj_info(idx,:);
        clearvars idx
    end

    
    % 3.3 calculate the in, out, tot strengths/degrees using the func countlinks in graphcompare
    if thr_h2==0 && norm==1
        
        [linksbase,~]=ins_countlinks(baseline,thr_h2);        % linkspre contains 1 row per channel, 1 column per window of h2
        OUTpre = median(linksbase.outstrength_norm,2);
        TOTpre = median(linksbase.totstrength_norm,2);
        %post1
        [linkspostA,~]=ins_countlinks(postA,thr_h2);
        OUTpostA = median(linkspostA.outstrength_norm,2);
        TOTpostA = median(linkspostA.totstrength_norm,2);
        %post2
        [linkspostB,~]=ins_countlinks(postB,thr_h2);
        OUTpostB = median(linkspostB.outstrength_norm,2);
        TOTpostB = median(linkspostB.totstrength_norm,2);
        %postAll
        [linkspostAll,~]=ins_countlinks(postAll,thr_h2);
        OUTpostAll = median(linkspostAll.outstrength_norm,2);
        TOTpostAll = median(linkspostAll.totstrength_norm,2);
        %sham
        [linkssham,~]=ins_countlinks(sham,thr_h2);
        OUTsham = median(linkssham.outstrength_norm,2);
        TOTsham = median(linkssham.totstrength_norm,2);

        if instim
            %stimA
            [linksstimA,~]=ins_countlinks(stimA,thr_h2);
            OUTstimA = median(linksstimA.outstrength_norm,2);
            TOTstimA = median(linksstimA.totstrength_norm,2);
    
            %stimB
            [linksstimB,~]=ins_countlinks(stimB,thr_h2);
            OUTstimB = median(linksstimB.outstrength_norm,2);
            TOTstimB = median(linksstimB.totstrength_norm,2);
        end
        
        for chan = 1: size(linksbase.totstrength,1)  
            clearvars statsa statsb statsham statspAll
            [~,~,statsham] = ranksum(linkssham.totstrength_norm(chan,:), linksbase.totstrength_norm(chan,:));
            ztotsham(chan) = statsham.zval;

            if strcmp(base, 'baseline')
                [~,~,statsa] = ranksum(linkspostA.totstrength_norm(chan,:), linksbase.totstrength_norm(chan,:)); %zvalue of norm is same of not norm (I also checked that)
                ztota(chan)  = statsa.zval;
                [~,~,statsb] = ranksum(linkspostB.totstrength_norm(chan,:), linksbase.totstrength_norm(chan,:));
                ztotb(chan)  = statsb.zval;
                [~,~,statspAll] = ranksum(linkspostAll.totstrength_norm(chan,:), linksbase.totstrength_norm(chan,:));
                ztotAll(chan)  = statspAll.zval;
                if instim
                    [~,~,statsastim] = ranksum(linksstimA.totstrength_norm(chan,:), linksbase.totstrength_norm(chan,:)); %zvalue of norm is same of not norm (I also checked that)
                    ztotstima(chan)  = statsastim.zval;
                    [~,~,statsbstim] = ranksum(linksstimB.totstrength_norm(chan,:), linksbase.totstrength_norm(chan,:));
                    ztotstimb(chan)  = statsbstim.zval;
                end
           
            elseif strcmp(base,'sham')
                [~,~,statsa] = ranksum(linkspostA.totstrength_norm(chan,:), linkssham.totstrength_norm(chan,:)); %zvalue of norm is same of not norm (I also checked that)
                ztota(chan)  = statsa.zval;
                [~,~,statsb] = ranksum(linkspostB.totstrength_norm(chan,:), linkssham.totstrength_norm(chan,:));
                ztotb(chan)  = statsb.zval;
                [~,~,statspAll] = ranksum(linkspostAll.totstrength_norm(chan,:), linkssham.totstrength_norm(chan,:));
                ztotAll(chan)  = statspAll.zval;
                if instim
                    [~,~,statsastim] = ranksum(linksstimA.totstrength_norm(chan,:), linkssham.totstrength_norm(chan,:)); %zvalue of norm is same of not norm (I also checked that)
                    ztotstima(chan)  = statsastim.zval;
                    [~,~,statsbstim] = ranksum(linksstimB.totstrength_norm(chan,:), linkssham.totstrength_norm(chan,:));
                    ztotstimb(chan)  = statsbstim.zval;
                end
                
            else
                error('Function only accepts either baseline or sham as base for zscore tests')
            end

        end
        
        TOTzA = ztota';
        TOTzB = ztotb';
        TOTzAll = ztotAll';
        TOTzSHAM = ztotsham';
        if instim
            TOTzstimA = ztotstima';
            TOTzstimB = ztotstimb';
        end

       clearvars ztota ztotb ztotsham ztotAll
        
  
        
    else %degrees
        [linksbase,~]=ins_countlinks(baseline,thr_h2);
        OUTpre = median(linksbase.outdegree,2);
        TOTpre = median(linksbase.totdegree,2);
        %post1
        [linkspostA,~]=ins_countlinks(postA,thr_h2);
        OUTpostA = median(linkspostA.outdegree,2);
        TOTpostA = median(linkspostA.totdegree,2);
        %post2
        [linkspostB,~]=ins_countlinks(postB,thr_h2);
        OUTpostB = median(linkspostB.outdegree,2);
        TOTpostB = median(linkspostB.totdegree,2);
        %postAll
        [linkspostAll,~]=ins_countlinks(postAll,thr_h2);
        OUTpostAll = median(linkspostAll.outdegree,2);
        TOTpostAll = median(linkspostAll.totdegree,2);
        %sham
        [linkssham,~]=ins_countlinks(sham,thr_h2);
        OUTsham = median(linkssham.outdegree,2);
        TOTsham = median(linkssham.totdegree,2);

        if instim
            %stimA
            [linksstimA,~]=ins_countlinks(stimA,thr_h2);
            OUTstimA = median(linksstimA.outdegree,2);
            TOTstimA = median(linksstimA.outdegree,2);
    
            %stimB
            [linksstimB,~]=ins_countlinks(stimB,thr_h2);
            OUTstimB = median(linksstimB.outdegree,2);
            TOTstimB = median(linksstimB.outdegree,2);
        end
        
        for chan = 1: size(linksbase.totstrength,1)  
            clearvars statsa statsb statsham statspAll
            [~,~,statsa] = ranksum(linksbase.totdegree(chan,:), linkspostA.totdegree(chan,:));
            ztota(chan) = statsa.zval;
            [~,~,statsb] = ranksum(linksbase.totdegree(chan,:), linkspostB.totdegree(chan,:));
            ztotb(chan) = statsb.zval;
            [~,~,statspAll] = ranksum(linksbase.totdegree(chan,:), linkspostAll.totdegree(chan,:));
            ztotAll(chan) = statspAll.zval;
            [~,~,statsham] = ranksum( linkssham.totdegree(chan,:), linksbase.totdegree(chan,:));
            ztotsham(chan) = statsham.zval;
            if instim
                [~,~,statsastim] = ranksum(linksstimA.totdegree(chan,:), linkssham.totdegree(chan,:)); %zvalue of norm is same of not norm (I also checked that)
                ztotstima(chan)  = statsastim.zval;
                [~,~,statsbstim] = ranksum(linksstimB.totdegree(chan,:), linkssham.totdegree(chan,:));
                ztotstimb(chan)  = statsbstim.zval;
            end
        end
        
        TOTzA = ztota';
        TOTzB = ztotb';
        TOTzAll = ztotAll';
        TOTzSHAM = ztotsham';
        if instim 
            TOTzstimA = ztotstima';
            TOTzstimB = ztotstimb';
        end

       clearvars ztota ztotb ztotsham statsbstim statsastim

    end
    
    if instim
        widedata = [widedata; string([subj_info OUTpre TOTpre OUTstimA TOTstimA OUTstimB TOTstimB OUTpostA TOTpostA OUTpostB TOTpostB OUTpostAll TOTpostAll OUTsham TOTsham TOTzstimA TOTzstimB TOTzA TOTzB TOTzAll TOTzSHAM])];
    else
        widedata = [widedata; string([subj_info OUTpre TOTpre OUTpostA TOTpostA OUTpostB TOTpostB OUTpostAll TOTpostAll OUTsham TOTsham TOTzA TOTzB TOTzAll TOTzSHAM])];
    end        
    
    clearvars info hvalues datapre datapostA datapostB subj channels channels_mono ...
        check linkspre linkspostA linkspostB subchan OUTpre TOTpre OUTpostA TOTpostA OUTpostB TOTpostB OUTsham TOTsham H subj_struc TOTzvalA TOTzvalB TOTzSHAM ...
        OUTstimA TOTstimA OUTstimB TOTstimB TOTzstimA ztotstima TOTzstimB ztotstimb TOTzAll OUTpostAll TOTpostAll
end

if instim
    FCtable = array2table(widedata, 'VariableNames',{varnames{:,:},'OUTpre' ,'TOTpre', ...
        'OUTstimA', 'TOTstimA', 'OUTstimB', 'TOTstimB', 'OUTpostA', 'TOTpostA', 'OUTpostB', ...
        'TOTpostB', 'OUTpostAll', 'TOTpostAll', 'OUTsham', 'TOTsham', 'TOTzstimA', 'TOTzstimB', 'TOTzA', 'TOTzB', 'TOTzAll', 'TOTzSHAM'});
else
    FCtable = array2table(widedata, 'VariableNames',{varnames{:,:},'OUTpre' ,...
        'TOTpre', 'OUTpostA', 'TOTpostA', 'OUTpostB', 'TOTpostB', 'OUTpostAll', ...
        'TOTpostAll', 'OUTsham', 'TOTsham', 'TOTzA', 'TOTzB', 'TOTzAll', 'TOTzSHAM'});
end   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Connectivity matrices figures etc
% format long g

% sorting by EF
if instim
    nsec = 7;
else 
    nsec = 5;
end

for index = 1:nsec:length(myfiles)
     
    baseline = load(myfiles(index).name); 
    postA = load(myfiles(index+1).name);
    postB = load(myfiles(index+2).name);
    sham = load(myfiles(index+3).name);
    postAll = load(myfiles(index+4).name);
    if instim
        stimA = load(myfiles(index+5).name);
        stimB = load(myfiles(index+6).name);
    end
    
    channels = string(baseline.electrode_names)';
    
    subj = string(extractBetween(myfiles(index).name, 'sub-', '_ses'));
    chaninfosubj = chan_info(strcmp(chan_info.subj, subj), :);
    
    channels_mono = extractBefore(string(channels), '-');
    
    for i=1:length(channels)
        if ~isempty(chaninfosubj.simEF(strcmp(channels_mono(i), chaninfosubj.channel)))
            channels(i,2) = chaninfosubj.simEF(strcmp(channels_mono(i), chaninfosubj.channel));
        else channels(i,2) = NaN;
        end
    end
    
    [~, idx] = sort(str2double(channels(:,2)), 'descend');
    channels_sorted = channels(idx,:);


% re-order .mat files by simEF
    
    if instim
         clearvars i
        data = {baseline, sham, postA, postB, postAll, stimA, stimB};
        for file =1:nsec
            data{file}.electrode_names = data{file}.electrode_names(idx);
            for i = 1: size(data{file}.aw_h2,3)
                data{file}.aw_h2(:,:,i) = data{file}.aw_h2(idx,idx,i);
                data{file}.aw_lag(:,:,i) = data{file}.aw_lag(idx,idx,i);
            end
        end
        
    else
        clearvars i
            data = {baseline, sham, postA, postB, postAll};
            for file =1:nsec
                data{file}.electrode_names = data{file}.electrode_names(idx);
                for i = 1: size(data{file}.aw_h2,3)
                    data{file}.aw_h2(:,:,i) = data{file}.aw_h2(idx,idx,i);
                    data{file}.aw_lag(:,:,i) = data{file}.aw_lag(idx,idx,i);
                end
            end
    end

    if sorted ==1
%         cd("C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\FC\h2\sorted_by_simEF_new")
        baseS  = data{1};
        shamS  = data{2};
        postAS = data{3}; 
        postBS = data{4};
        postAllS = data{5};
        if instim
            stimAS = data{6};
            stimBS = data{7};
        end
%         save(strcat(subj,'_base_broad.mat'), '-struct', 'baseS')
%         save(strcat(subj,'_sham_broad.mat'),  '-struct','shamS')
%         save(strcat(subj,'_postA_broad.mat'), '-struct', 'postAS')
%         save(strcat(subj,'_postB_broad.mat'), '-struct', 'postBS')
%         if instim
%             save(strcat(subj,'_stimA_broad.mat'), '-struct', 'stimAS')
%             save(strcat(subj,'_stimB_broad.mat'), '-struct', 'stimBS')
%         end
    
%         cd(dir_data)
    end


%% GRAPHS
    if graph == 1
    
    % CALCULATE zscore FOR EACH ELEMENT of conn matrix across z direction 
    % first line: post vs baseline, second line: post vs sham
    if sorted
        for i = 1:length(channels)
        for j = 1:length(channels)
        [~,~,zvalsham] = ranksum(squeeze(shamS.aw_h2(i, j, :)), squeeze(baseS.aw_h2(i, j, :)));
        zsham(i, j) = zvalsham.zval;

        [~,~,zvalpostA] = ranksum(squeeze(postAS.aw_h2(i, j, :)), squeeze(baseS.aw_h2(i, j, :)));
        zpostA(i, j) = zvalpostA.zval;
        [~,~,zvalpostB] = ranksum(squeeze(postBS.aw_h2(i, j, :)), squeeze(baseS.aw_h2(i, j, :)));
        zpostB(i, j) = zvalpostB.zval;
        [~,~,zvalpostAll] = ranksum(squeeze(postAllS.aw_h2(i, j, :)), squeeze(baseS.aw_h2(i, j, :)));
        zpostAll(i, j) = zvalpostAll.zval;

        [~,~,zvalpostAsham] = ranksum(squeeze(postAS.aw_h2(i, j, :)), squeeze(shamS.aw_h2(i, j, :)));
        zpostAsham(i, j) = zvalpostAsham.zval;
        [~,~,zvalpostBsham] = ranksum(squeeze(postBS.aw_h2(i, j, :)), squeeze(shamS.aw_h2(i, j, :)));
        zpostBsham(i, j) = zvalpostBsham.zval;
        [~,~,zvalpostAllsham] = ranksum(squeeze(postAllS.aw_h2(i, j, :)), squeeze(shamS.aw_h2(i, j, :)));
        zpostAllsham(i, j) = zvalpostAllsham.zval;

        if instim
            [~,~,zvalstimA] = ranksum(squeeze(stimAS.aw_h2(i, j, :)), squeeze(baseS.aw_h2(i, j, :)));
            zstimA(i, j) = zvalstimA.zval;
            [~,~,zvalstimB] = ranksum(squeeze(stimBS.aw_h2(i, j, :)), squeeze(baseS.aw_h2(i, j, :)));
            zstimB(i, j) = zvalstimB.zval;

            [~,~,zvalstimAsham] = ranksum(squeeze(stimAS.aw_h2(i, j, :)), squeeze(shamS.aw_h2(i, j, :)));
            zstimAsham(i, j) = zvalstimAsham.zval;
            [~,~,zvalstimBsham] = ranksum(squeeze(stimBS.aw_h2(i, j, :)), squeeze(shamS.aw_h2(i, j, :)));
            zstimBsham(i, j) = zvalstimBsham.zval;
        end
        end
        end

    else

        for i = 1:length(channels)
        for j = 1:length(channels)
        [~,~,zvalsham] = ranksum(squeeze(sham.aw_h2(i, j, :)), squeeze(baseline.aw_h2(i, j, :)));
        zsham(i, j) = zvalsham.zval;
        [~,~,zvalpostA] = ranksum(squeeze(postA.aw_h2(i, j, :)), squeeze(baseline.aw_h2(i, j, :)));
        zpostA(i, j) = zvalpostA.zval;
        [~,~,zvalpostB] = ranksum(squeeze(postB.aw_h2(i, j, :)), squeeze(baseline.aw_h2(i, j, :)));
        zpostB(i, j) = zvalpostB.zval;
        [~,~,zvalpostAll] = ranksum(squeeze(postAll.aw_h2(i, j, :)), squeeze(baseline.aw_h2(i, j, :)));
        zpostAll(i, j) = zvalpostAll.zval;

        [~,~,zvalpostAsham] = ranksum(squeeze(postA.aw_h2(i, j, :)), squeeze(sham.aw_h2(i, j, :)));
        zpostAsham(i, j) = zvalpostAsham.zval;
        [~,~,zvalpostBsham] = ranksum(squeeze(postB.aw_h2(i, j, :)), squeeze(sham.aw_h2(i, j, :)));
        zpostBsham(i, j) = zvalpostBsham.zval;
        [~,~,zvalpostAllsham] = ranksum(squeeze(postAll.aw_h2(i, j, :)), squeeze(sham.aw_h2(i, j, :)));
        zpostAllsham(i, j) = zvalpostAllsham.zval;

        if instim
            [~,~,zvalstimA] = ranksum(squeeze(stimA.aw_h2(i, j, :)), squeeze(baseline.aw_h2(i, j, :)));
            zstimA(i, j) = zvalstimA.zval;
            [~,~,zvalstimB] = ranksum(squeeze(stimB.aw_h2(i, j, :)), squeeze(baseline.aw_h2(i, j, :)));
            zstimB(i, j) = zvalstimB.zval;
            [~,~,zvalstimAsham] = ranksum(squeeze(stimA.aw_h2(i, j, :)), squeeze(sham.aw_h2(i, j, :)));
            zstimAsham(i, j) = zvalstimAsham.zval;
            [~,~,zvalstimBsham] = ranksum(squeeze(stimB.aw_h2(i, j, :)), squeeze(sham.aw_h2(i, j, :)));
            zstimBsham(i, j) = zvalstimBsham.zval;
        end
        end
        end
    end

    
    
        base_mat  = mean(baseline.aw_h2,3);
        sham_mat  = mean(sham.aw_h2,3);
        postA_mat = mean(postA.aw_h2,3);
        postB_mat = mean(postB.aw_h2,3);
        postAll_mat = mean(postAll.aw_h2,3);
        
        base_sort  = base_mat(idx,idx);
        sham_sort  = sham_mat(idx,idx);
        postA_sort = postA_mat(idx,idx);
        postB_sort = postB_mat(idx,idx);
        postAll_sort = postAll_mat(idx,idx);
        
       

        channels_sorted(ismissing(channels_sorted(:,2)),2) = "missing";
        simEFvals = extractBefore(channels_sorted(:,2), 5);
        

        %% h2 conn matrices
        lim = [-4 4]; 
if sorted
        nrow =2; 
        if instim
            ncol = 5; nsubf = 0;
        else 
            ncol = 3; nsubf = -2;
        end

        color = viridis;
% control of sham vs baseline for each subjects

        figure('Name', strcat('Z-values of FC ordered by simEF || sub-', subj));
        imagesc(zsham)
        colormap(color)
        title("zvalue sham vs baseline")
        xticks(1:length(channels));
        xticklabels(strcat(channels_sorted(:,1), '/', simEFvals));
        xtickangle(60)
        yticks(1:length(channels));
        yticklabels(channels_sorted(:,1));
        ytickangle(0)
        colorbar
        caxis(lim);
        ylabel = "channel";
        xlabel = "channel";

% 1 figure per subject with postA and B (+stimA/B if instim). 1st row: ranksum vs baseline, 2nd row: ranksum vs sham

        figure('Name', strcat('Z-values of FC ordered by simEF || sub-', subj));

        if instim      
            subplot(nrow, ncol,1)
            imagesc(zstimA)
            colormap(color)
            title("zvalue stimA vs baseline")
            xticks(1:length(channels));
            xticklabels(strcat(channels_sorted(:,1), '/', simEFvals));
            xtickangle(60)
            yticks(1:length(channels));
            yticklabels(channels_sorted(:,1));
            ytickangle(0)
            colorbar
            caxis(lim);
            ylabel = "channel";
            xlabel = "channel";
    
            subplot(nrow, ncol, 2)
            imagesc(zstimB)
            colormap(color)
            title("zvalue stimB vs baseline")
            xticks(1:length(channels));
            xticklabels(strcat(channels_sorted(:,1), '/', simEFvals));
            xtickangle(60)
            yticks(1:length(channels));
            yticklabels(channels_sorted(:,1));
            ytickangle(0)
            colorbar
            caxis(lim);
            ylabel = "channel";
            xlabel = "channel";

        end

        subplot(nrow, ncol,3 + nsubf)
        imagesc(zpostA)
        colormap(color)
        title("zvalue postA vs baseline")
        xticks(1:length(channels));
        xticklabels(strcat(channels_sorted(:,1), '/', simEFvals));
        xtickangle(60)
        yticks(1:length(channels));
        yticklabels(channels_sorted(:,1));
        ytickangle(0)
        colorbar
        caxis(lim);
        ylabel = "channel";
        xlabel = "channel";
        
        subplot(nrow, ncol,4  + nsubf)
        imagesc(zpostB)
        colormap(color)
        title("zvalue postB vs baseline")
        xticks(1:length(channels));
        xticklabels(strcat(channels_sorted(:,1), '/', simEFvals));
        xtickangle(60)
        yticks(1:length(channels));
        yticklabels(channels_sorted(:,1));
        ytickangle(0)
        colorbar
        caxis(lim);
        ylabel = "channel";
        xlabel = "channel";
        
        subplot(nrow, ncol, 5 + nsubf)
        imagesc(zpostAll)
        colormap(color)
        title("zvalue postAll vs baseline")
        xticks(1:length(channels));
        xticklabels(strcat(channels_sorted(:,1), '/', simEFvals));
        xtickangle(60)
        yticks(1:length(channels));
        yticklabels(channels_sorted(:,1));
        ytickangle(0)
        colorbar
        caxis(lim);
        ylabel = "channel";
        xlabel = "channel";
                
% vs sham

    if instim
        subplot(nrow, ncol, 6)
            imagesc(zstimAsham)
            colormap(color)
            title("zvalue stimA vs sham")
            xticks(1:length(channels));
            xticklabels(strcat(channels_sorted(:,1), '/', simEFvals));
            xtickangle(60)
            yticks(1:length(channels));
            yticklabels(channels_sorted(:,1));
            ytickangle(0)
            colorbar
            caxis(lim);
            ylabel = "channel";
            xlabel = "channel";
            
            subplot(nrow, ncol, 7)
            imagesc(zstimBsham)
            colormap(color)
            title("zvalue stimB vs sham")
            xticks(1:length(channels));
            xticklabels(strcat(channels_sorted(:,1), '/', simEFvals));
            xtickangle(60)
            yticks(1:length(channels));
            yticklabels(channels_sorted(:,1));
            ytickangle(0)
            colorbar
            caxis(lim);
            ylabel = "channel";
            xlabel = "channel";
    end

        subplot(nrow, ncol, 8 + nsubf*2)
        imagesc(zpostAsham)
        colormap(color)
        title("zvalue postA vs sham")
        xticks(1:length(channels));
        xticklabels(strcat(channels_sorted(:,1), '/', simEFvals));
        xtickangle(60)
        yticks(1:length(channels));
        yticklabels(channels_sorted(:,1));
        ytickangle(0)
        colorbar
        caxis(lim);
        ylabel = "channel";
        xlabel = "channel";
        
        subplot(nrow, ncol, 9 + nsubf*2)
        imagesc(zpostBsham)
        colormap(color)
        title("zvalue postB vs sham")
        xticks(1:length(channels));
        xticklabels(strcat(channels_sorted(:,1), '/', simEFvals));
        xtickangle(60)
        yticks(1:length(channels));
        yticklabels(channels_sorted(:,1));
        ytickangle(0)
        colorbar
        caxis(lim);
        ylabel = "channel";
        xlabel = "channel";

        subplot(nrow, ncol, 10 + nsubf*2)
        imagesc(zpostAllsham)
        colormap(color)
        title("zvalue postAll vs sham")
        xticks(1:length(channels));
        xticklabels(strcat(channels_sorted(:,1), '/', simEFvals));
        xtickangle(60)
        yticks(1:length(channels));
        yticklabels(channels_sorted(:,1));
        ytickangle(0)
        colorbar
        caxis(lim);
        ylabel = "channel";
        xlabel = "channel";


        clearvars channels_sorted channels zsham zstimA zstimB zpostA zpostB zpostAll zstimAsham zstimBsham zpostAsham zpostBsham zpostAllsham...
            sham postA postB postAll stimA stimB shamS postAS postBS postAllS stimAS stimBS baseline baseS...
            zvalsham zvalpostA zvalpostB zvalstimB zvalstimA zvalpostAll zvalpostAsham zvalpostBsham zvaalstimAsham zvalstimBsham zvalpostAllsham
end
%         %% LAG
%         base_lag = mean(baseline.aw_lag,3);
%         sham_lag = mean(sham.aw_lag,3);
%         postA_lag = mean(postA.aw_lag,3);
%         postB_lag = mean(postB.aw_lag,3);
%         
%         baseL_sort = base_lag(idx,idx);
%         shamL_sort = sham_lag(idx,idx);
%         postAL_sort = postA_lag(idx,idx);
%         postBL_sort = postB_lag(idx,idx);
%         
%         %% LAG conn matrices
%         
%         lim = [-22 22];
%         map = "viridis";
%         
%         figure('Name', subj);
%         subplot(2,2,1)
%         imagesc(baseL_sort)
%         colormap(map)
%         title("baseline")
%         xticks(1:length(channels));
%         xticklabels(channels_sorted(:,1));
%         xtickangle(60)
%         yticks(1:length(channels));
%         yticklabels(channels_sorted(:,1));
%         ytickangle(0)
%         colorbar
%         caxis(lim);
%         ylabel = "channel";
%         xlabel = "channel";
%         
%         subplot(2,2,2)
%         imagesc(shamL_sort)
%         colormap(map)
%         title("sham")
%         xticks(1:length(channels));
%         xticklabels(channels_sorted(:,1));
%         xtickangle(60)
%         yticks(1:length(channels));
%         yticklabels(channels_sorted(:,1));
%         ytickangle(0)
%         colorbar
%         caxis(lim);
%         ylabel = "channel";
%         xlabel = "channel";
%         
%         subplot(2,2,3)
%         imagesc(postAL_sort)
%         colormap(map)
%         title("postA")
%         xticks(1:length(channels));
%         xticklabels(channels_sorted(:,1));
%         xtickangle(60)
%         yticks(1:length(channels));
%         yticklabels(channels_sorted(:,1));
%         ytickangle(0)
%         colorbar
%         caxis(lim);
%         ylabel = "channel";
%         xlabel = "channel";
%         
%         subplot(2,2,4)
%         imagesc(postBL_sort)
%         colormap(map)
%         title("postB")
%         xticks(1:length(channels));
%         xticklabels(channels_sorted(:,1));
%         xtickangle(60)
%         yticks(1:length(channels));
%         yticklabels(channels_sorted(:,1));
%         ytickangle(0)
%         colorbar
%         caxis(lim);
%         ylabel = "channel";
%         xlabel = "channel";

        clearvars channels_sorted channels zsham zpostA zpostAsham zpostBsham zvalpostBsham zvalpostAsham zvalpostB zvalpostA zvalsham

        
    end

end

end


