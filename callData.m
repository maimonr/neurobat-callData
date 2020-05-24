classdef callData
    % CallData Version 2.0 Maimon Rose 06/30/17
    % Includes the following experiment types: Maimon's ephys, Deafened bats,
    % Lesioned bats, and Tobias' automatic training.
    %
    % Construct a callData object as follows:
    %
    % cData = callData(fs,expType)
    % fs - audio sampling rate in Hz
    % expType - string correspoding to experiment type
    
    properties(SetAccess = public)
        baseDirs % directories to search for call files
        batNums % bat ID numbers (as strings)
        dateFormat % format for importing dates as datetimes
        birthDates % DOBs for all bats
        treatment
        treatmentType
        nBats % number of bats
        cepstralFlag = true % perform cepstral analysis?
        spCorrFlag = false % perform spCorr analysis?
        batName %name of bat pairs in training
        batchNum % for experiments in separate batches
        callNum %call number (subset of recNum)
        recNum %recording number for the training session
        callTrigger %dt for autTrain
        callType %'callOnly' or 'callTrigger'
        sessionType %experiment type in the autoTrain
        sessionID %numerical counter for each sessionTye
        micType %which box/microphone used
        xrun %dropped sample number
        callID %unique ID for each call
        callEcho % flag to load either calls or echolocations
        loggerNum % audio logger serial number
        manual_call_class % variable to store potenial manual classification
        
        fs % sampling rate
        expType % String indicating which experiment we are dealing with
        exp_session_type % String indicating which session of an experiment we are dealing with
        
        callWF % call waveform
        expDay % experimental date (in datetime format)
        exp_date_num % experimental date (in datenum format, for increase computation speed)
        recTime % time of recording (in datetime format)
        daysOld % number of days from birth for this call
        batNum % bat ID number correspoding to this call
        fName % file path to file from which call was cut
        fName_cut % file path to file with cut call
        nCalls % total number of calls
        callLength % length of call
        file_call_pos % % (for 'ephys') position in samples during session
        callPos % (for 'ephys') position in samples during session
        yinF0 % 'best' fundamental freq. calculated using the yin algorithm
        spCorrF0 % mean fundamental freq. calculated using the spCorr algorithm
        weinerEntropy % mean WE
        spectralEntropy % mean SE
        centroid % mean SE
        energyEntropy % mean EE
        RMS % mean RMS
        ap0 % 'best' coarse aperiodicity
        pitchGoodness
        cepstralF0 % fundamental freq. calculated using a cepstral algorithm
        
        % parameters as above, broken down into windows during the call:
        yinF0Win
        ap0Win
        weinerEntropyWin
        spectralEntropyWin
        centroidWin
        energyEntropyWin
        RMSWin
        pitchGoodnessWin
        cepstralF0Win
        spCorrF0Win
        
        maxCalls = 150e3
        minF0 = 200 % in Hz
        maxF0 = 20e3 % in Hz
        
        windowLength = 1.5e-3 % in ms
        integrationWindow = 1.5e-3
        overlap = 1e-3 % in ms
        thresh = 0.01
        yinF0_interp = false
        EW_mic_HP_freq = 200
        tapers % precalculate tapers for mt spectral estimation
        
        loadWF = true % Load in all call waveforms or not
        compensation % compensation impulse response filter
    end
    properties (Dependent)
        oldestAge
        youngestAge
        yinParams
        featureParams
    end
    methods
        function cData = callData(fs, varargin) % Function to initialize, load, and perform calculations
            cData.fs = fs;
            
            pnames = {'expType','exp_session_type','callEcho'};
            dflts  = {'adult','communication','call'};
            [expType,exp_session_type,callEcho] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            
            cData.exp_session_type = exp_session_type;
            cData.expType = expType;
            cData.callEcho = callEcho;
            
            % initialize all parameters:
            
            cData = load_call_data(cData);
            
            [cData.yinF0, cData.spCorrF0,...
                cData.weinerEntropy, cData.spectralEntropy, cData.centroid,...
                cData.energyEntropy, cData.RMS, cData.ap0,cData.pitchGoodness,...
                cData.cepstralF0] = deal(zeros(cData.nCalls,1));
            
            [cData.yinF0Win, cData.ap0Win, cData.weinerEntropyWin,...
                cData.spectralEntropyWin, cData.centroidWin,cData.energyEntropyWin,...
                cData.RMSWin, cData.pitchGoodnessWin,cData.cepstralF0Win]...
                = deal(cell(cData.nCalls,1));
            
            cData.tapers = dpss(cData.windowLength*cData.fs,1.5);
            
            if ~cData.loadWF
                cData.callLength = zeros(cData.nCalls,1);
            end
            
            % iterate through all calls
            lastProgress = 0;
            
            for call_k = 1:cData.nCalls
                
                if cData.loadWF % if we already loaded all the data into this object
                    callWF = cData.callWF{call_k};
                else % if we need to load this data on the fly
                    if ~isempty(cData.fName{call_k})
                        callWF = loadCallWF_onTheFly(cData,call_k);
                        cData.callLength(call_k) = length(callWF)/cData.fs;
                    end
                end
                
                % Calculate windowed call features
                
                [f0,ap0,WE,SE,Cent,EE,Rms,PG,CepF0,spF0] = calculate_call_features(cData,callWF);
                
                cData.yinF0Win{call_k} = f0;
                cData.ap0Win{call_k} = ap0;
                
                % store those features
                
                [cData.weinerEntropyWin{call_k},cData.spectralEntropyWin{call_k},...
                    cData.centroidWin{call_k},cData.energyEntropyWin{call_k},...
                    cData.RMSWin{call_k},cData.pitchGoodnessWin{call_k},...
                    cData.cepstralF0Win{call_k}, cData.spCorrF0Win{call_k}]...
                    = deal(WE,SE,Cent,EE,Rms,PG,CepF0,spF0);
                
                % store average values of those features
                
                cData.yinF0(call_k) = nanmedian(f0);
                cData.ap0(call_k) = nanmedian(ap0);
                
                cData.weinerEntropy(call_k) = mean(cData.weinerEntropyWin{call_k});
                cData.spectralEntropy(call_k) = mean(cData.spectralEntropyWin{call_k});
                cData.centroid(call_k) = mean(cData.centroidWin{call_k});
                cData.energyEntropy(call_k) = mean(cData.energyEntropyWin{call_k});
                cData.RMS(call_k) = mean(cData.RMSWin{call_k});
                
                if cData.cepstralFlag % if we want to use the cepstral algorithm
                    cData.pitchGoodness(call_k) = nanmedian(cData.pitchGoodnessWin{call_k});
                    cData.cepstralF0(call_k) = nanmedian(cData.cepstralF0Win{call_k});
                end
                
                if cData.spCorrFlag
                    cData.spCorrF0(call_k) = nanmedian(cData.spCorrF0(call_k));
                end
                
                progress = 100*(call_k/cData.nCalls);
                
                if mod(progress,10) < mod(lastProgress,10)
                    fprintf('%d %% of calls processed\n',round(progress));
                end
                
                lastProgress = progress;
                
            end
            
            cData.exp_date_num = datenum(cData.expDay);
            
        end
        function [f0,ap0,WE,SE,Cent,EE,Rms,PG,CepF0,spF0] = calculate_call_features(cData,callWF)
            
            nSamples = length(callWF);
            YP = cData.yinParams;
            YP.bufsize = min(nSamples,cData.fs*2);
            YP.range = [1 nSamples];
            [f0, ap0] = calculate_yin_F0(callWF,YP,cData.yinF0_interp);
            
            [WE,SE,Cent,EE,Rms,PG,CepF0,spF0] = getCallFeatures(callWF,cData.featureParams);
            
        end
        function n = numArgumentsFromSubscript(obj,~,~)
            n = numel(obj);
        end
        function cData = load_call_data(cData,varargin)
            
            if any(strcmp(cData.expType,{'juvenile','adult','adult_operant'})) % Maimon's ephys experiments
                
                eData = ephysData(cData.expType);
                cData.loadWF = false;
                
                irc_fname = fullfile('E:\ephys\avisoft_mic1_calibration','Comp_IR_avisoft1_freefield_20161109T100619.mat');
                
                if ~strcmp(cData.exp_session_type,'operant') && exist(irc_fname,'file')
                    irc = load(irc_fname);
                    cData.compensation.irc = irc.irc;
                    cData.compensation.fs = cData.fs;
                    cData.fs = irc.fs;
                else
                    disp('No compensation filter found')
                end
                
                cData.baseDirs = eData.baseDirs;
                cData.batNums = eData.batNums;
                cData.dateFormat = eData.dateFormat;
                cData.birthDates = eData.birthDates;
                cData.nBats = length(cData.batNums);
                
                [cData.callWF, cData.batNum, cData.fName] = deal(cell(cData.maxCalls,1));
                [cData.daysOld, cData.callLength] = deal(zeros(cData.maxCalls,1));
                cData.callPos = zeros(cData.maxCalls,2);
                cData.file_call_pos = zeros(cData.maxCalls,2);
                cData.expDay = datetime([],[],[]);
                
                call_k = 1;
                
                all_cut_call_data = get_cut_call_data(cData);
                
                for k = 1:length(all_cut_call_data)
                    cut_call_data = all_cut_call_data{k};
                    if ~isempty(cut_call_data) && ~isempty({cut_call_data.fName})
                        
                        cutCalls = {cut_call_data.cut};
                        call_pos_expDay = vertcat(cut_call_data.corrected_callpos)/1e3; % convert to seconds
                        file_call_pos_expDay = vertcat(cut_call_data.callpos);
                        call_length_expDay = cellfun(@length, cutCalls)/cData.fs;
                        
                        if any(isnat([cut_call_data.expDay]))
                            try
                                expDays = arrayfun(@(x) datetime(strrep(regexp(x.fName,'\\\d{8}\','match'),'\',''),'inputFormat','yyyyMMdd'),cut_call_data,'un',0);
                            catch
                                expDays = arrayfun(@(x) datetime(strrep(regexp(x.fName,'\\\d{8}\','match'),'\',''),'inputFormat','MMddyyyy'),cut_call_data,'un',0);
                            end
                            [cut_call_data(:).expDay] = deal(expDays{:});
                        end
                        
                        assert(all(isdatetime([cut_call_data.expDay])) && length(unique([cut_call_data.expDay])) == 1)
                        expDatetime = cut_call_data(1).expDay;
                        
                        for c = 1:length(cutCalls)
                            if ~cut_call_data(c).noise && ~any(isnan(cut_call_data(c).corrected_callpos))
                                cData.callLength(call_k) = call_length_expDay(c);
                                cData.callPos(call_k,:) = call_pos_expDay(c,:);
                                cData.file_call_pos(call_k,:) = file_call_pos_expDay(c,:);
                                cData.expDay(call_k,1) = expDatetime;
                                cData.batNum{call_k} = cut_call_data(c).batNum;
                                if length(cData.batNum{call_k}) == 1
                                    b = strcmp(cData.batNums,cData.batNum{call_k});
                                    if any(b)
                                        cData.daysOld(call_k) = days(expDatetime - cData.birthDates{b});
                                    else
                                        cData.daysOld(call_k) = NaN;
                                    end
                                end
                                current_fname = cut_call_data(c).fName;
                                if (ischar(current_fname) && exist(cut_call_data(c).fName,'file')) || (iscell(current_fname) && any(logical(cellfun(@exist,current_fname))))
                                    cData.fName{call_k} = current_fname;
                                else
                                    cData.fName{call_k} = strrep(cut_call_data(c).fName,'Z:\users\Maimon\ephys\','E:\ephys\juvenile_recording\');
                                    if ~logical(exist(cData.fName{call_k} ,'file'))
                                        continue
                                    end
                                end
                                cData.callID(call_k) = cut_call_data(c).uniqueID;
                                call_k = call_k + 1;
                            end
                        end
                    end
                end
                
                cData.nCalls = call_k-1;
                
                
            elseif strcmp(cData.expType,'lesion')
                
                maxBatchNums = 5;
                batchNums = 1:maxBatchNums;
                cData.loadWF = false;
                s = load('E:\lesion_recordings\all_lesion_bat_info.mat','all_lesion_bat_info');
                batInfo = s.all_lesion_bat_info;
                
                cData.dateFormat = 'yyyyMMdd';
                dateRegExp = '_\d{8}T';
                
                [cData.callWF, cData.fName, cData.treatment] = deal(cell(cData.maxCalls,1));
                [cData.callLength, cData.daysOld, cData.batchNum] = deal(zeros(cData.maxCalls,1));
                cData.expDay = datetime([],[],[]);
                cData.callID = [1:cData.maxCalls]';
                
                treatmentGroups = fieldnames(batInfo)';
                
                if ~isempty(varargin)
                    varargin = varargin{1};
                    for subs_k = 1:2:length(varargin)
                        switch varargin{subs_k}
                            
                            case 'treatment'
                                treatmentGroups = treatmentGroups(ismember(treatmentGroups,varargin{subs_k+1}));
                            case 'batchNum'
                                batchNums = varargin{subs_k+1};
                        end
                        
                    end
                end
                
                
                call_k = 1;
                for t = 1:length(treatmentGroups)
                    groupStr = treatmentGroups{t};
                    batchNumsGroup = find(~structfun(@isempty,batInfo.(groupStr).birthDates))';
                    batchNumsGroup = batchNumsGroup(ismember(batchNumsGroup,batchNums));
                    
                    for b = batchNumsGroup
                        bStr = ['batch' num2str(b)];
                        avg_birth_date = mean(batInfo.(groupStr).birthDates.(bStr));
                        callFiles = dir([batInfo.(groupStr).baseDir  batInfo.(groupStr).treatmentType.(bStr) '*' bStr '*_Call_*.mat']);
                        if ~isempty(callFiles)
                            for c = 1:length(callFiles)
                                cData.treatment{call_k} = batInfo.(groupStr).treatmentType.(bStr);
                                exp_date_str = regexp(callFiles(c).name,dateRegExp,'match');
                                exp_date_str = exp_date_str{1}(2:end-1);
                                cData.expDay(call_k) = datetime(exp_date_str,'inputFormat',cData.dateFormat);
                                if isdatetime(avg_birth_date)
                                    cData.daysOld(call_k) = days(cData.expDay(call_k) - avg_birth_date);
                                end
                                cData.fName{call_k} = [batInfo.(groupStr).baseDir callFiles(c).name];
                                cData.batchNum(call_k) = b;
                                call_k = call_k + 1;
                            end
                        end
                    end
                end
                cData.nCalls = call_k-1;
                
                
            elseif strcmp(cData.expType,'deafened')
                maxBatchNums = 5;
                batchNums = 1:maxBatchNums;
                cData.loadWF = false;
                s = load('E:\deafened_recordings\all_deaf_bat_info.mat','all_deaf_bat_info');
                batInfo = s.all_deaf_bat_info;
                
                cData.dateFormat = 'yyyyMMddHHmmss';
                dateRegExp = '_\d{8}T';
                timeRegExp = 'T\d{6}_';
                
                [cData.callWF, cData.fName, cData.fName_cut, cData.treatment] = deal(cell(cData.maxCalls,1));
                [cData.callLength, cData.daysOld, cData.batchNum] = deal(zeros(cData.maxCalls,1));
                cData.expDay = datetime([],[],[]);
                cData.callID = [1:cData.maxCalls]';
                cData.file_call_pos = zeros(cData.maxCalls,2);
                
                treatmentGroups = fieldnames(batInfo)';
                
                if ~isempty(varargin)
                    varargin = varargin{1};
                    for subs_k = 1:2:length(varargin)
                        switch varargin{subs_k}
                            
                            case 'treatment'
                                treatmentGroups = treatmentGroups(ismember(treatmentGroups,varargin{subs_k+1}));
                            case 'batchNum'
                                batchNums = varargin{subs_k+1};
                        end
                        
                    end
                end
                
                
                call_k = 1;
                for t = 1:length(treatmentGroups)
                    groupStr = treatmentGroups{t};
                    batchNumsGroup = find(~structfun(@isempty,batInfo.(groupStr).birthDates));
                    batchNumsGroup = batchNumsGroup(ismember(batchNumsGroup,batchNums))';
                    group_data_baseDir = batInfo.(groupStr).baseDir(1:strfind(batInfo.deaf.baseDir,'all_cut_calls')-1);
                    
                    for b = batchNumsGroup
                        bStr = ['batch' num2str(b)];
                        avg_birth_date = mean(batInfo.(groupStr).birthDates.(bStr));
                        callFiles = dir([batInfo.(groupStr).baseDir  batInfo.(groupStr).treatmentType.(bStr) '*' bStr '*_Call_*.mat']);
                        raw_rec_dir = [group_data_baseDir 'all_raw_recordings\' batInfo.(groupStr).data_dir_str.(bStr) filesep];
                        if ~isempty(callFiles)
                            for c = 1:length(callFiles)
                                cData.treatment{call_k} = batInfo.(groupStr).treatmentType.(bStr);
                                exp_date_str = regexp(callFiles(c).name,dateRegExp,'match');
                                time_date_str = regexp(callFiles(c).name,timeRegExp,'match');
                                exp_datetime_str = [exp_date_str{1}(2:end-1) time_date_str{1}(2:end-1)];
                                cData.expDay(call_k) = datetime(exp_datetime_str,'inputFormat',cData.dateFormat);
                                if isdatetime(avg_birth_date)
                                    cData.daysOld(call_k) = days(cData.expDay(call_k) - avg_birth_date);
                                end
                                cData.fName_cut{call_k} = [batInfo.(groupStr).baseDir callFiles(c).name];
                                cData.fName{call_k} = [raw_rec_dir callFiles(c).name(1:strfind(callFiles(c).name,'_Call_')-1) '.mat'];
                                if ~exist(cData.fName{call_k},'file')
                                    cData.fName{call_k} = [];
                                end
                                cData.batchNum(call_k) = b;
                                s = load(cData.fName{call_k},'callpos');
                                cData.file_call_pos(call_k,:) = s.callpos;
                                call_k = call_k + 1;
                            end
                        end
                    end
                end
                cData.nCalls = call_k-1;
                cData.expDay = cData.expDay';
            elseif strcmp(cData.expType,'autoTrain')
                cData.loadWF = false;
                cData.baseDirs = 'C:\Users\tobias\Desktop\analysis\bataudio\call groups\all\';
                %cData.baseDirs = 'C:\Users\tobias\Desktop\analysis\bataudio\April2017\cut_preprocessed\';
                %cData.baseDirs = 'C:\Users\tobias\Desktop\analysis\bataudio\test\';
                callFiles = dir([cData.baseDirs '*.mat']);
                cData.nCalls = length(callFiles);
                [cData.fName, cData.batName, cData.callType, cData.micType, cData.sessionType] = deal(cell(cData.nCalls,1)); % initialize call data cells
                [cData.recNum, cData.callNum, cData.xrun, cData.sessionID,] = deal(zeros(cData.nCalls,1)); % initialize call data arrays
                cData.expDay = datetime([],[],[]);
                for call_k = 1:length(callFiles)
                    s = load([cData.baseDirs callFiles(call_k).name]);
                    cData.fName{call_k} = callFiles(call_k).name;
                    cData.batName{call_k} = s.batName;
                    cData.sessionType{call_k} = s.sessionType;
                    %cData.callWF{call_k} = s.rawData';
                    cData.callNum(call_k) = s.callNum;
                    cData.recNum(call_k) = s.recNum;
                    sessID = s.sessionID;
                    if ischar(sessID)
                        sessID = str2double(sessID);
                    end
                    if size(sessID) == [1 2]
                        sessID = 0;
                    end
                    cData.sessionID(call_k) = sessID;
                    if isfield(s,'xrun')
                        cData.xrun(call_k) = s.xrun;
                    else
                        cData.xrun(call_k) = 0;
                    end
                    if isfield(s,'callTrigger')
                        cData.expDay(call_k) = s.callTrigger;
                    else
                        cData.expDay(call_k) = s.dt;
                    end
                    cData.callType{call_k} = s.callType;
                    cData.micType{call_k} = s.micType;
                end
                cData.expDay = cData.expDay';
                
            elseif strcmp(cData.expType,'pratData')
                cData.loadWF = false;
                cData.baseDirs = 'E:\Yossi_vocalization_data\';
                
                [batInfo,batMetadata] = get_yossi_bat_info(cData.baseDirs);
                
                treatmentTypes = {'isolated','control'};
                
                cData.nCalls = sum(cellfun(@(x) length(vertcat(batInfo.(x)(:).callPos)),treatmentTypes));
                [cData.fName, cData.treatment, cData.batNum] = deal(cell(cData.nCalls,1));
                [cData.callLength, cData.daysOld, cData.treatment] = deal(zeros(cData.nCalls,1));
                cData.recTime = datetime(zeros(0,3));
                cData.batNums = cellfun(@(x) unique({batInfo.(x).batNum}),treatmentTypes,'un',0);
                cData.batNums = [cData.batNums{:}];
                cData.file_call_pos = zeros(cData.nCalls,2);
                
                cData.birthDates = datetime(zeros(0,3));
                for b = 1:length(cData.batNums)
                    dob = batMetadata.DOB{batMetadata.ID==str2double(cData.batNums{b})};
                    cData.birthDates(b) = dob;
                end
                cData.callID = [1:cData.nCalls]';
                call_k = 1;
                for t = treatmentTypes
                    treat = t{1};
                    for file_k = 1:length(batInfo.(treat))
                        for file_call_k = 1:size(batInfo.(treat)(file_k).callPos,1)
                            cData.treatment{call_k} = treat;
                            cData.fName{call_k} = [cData.baseDirs batInfo.(treat)(file_k).fName];
                            cData.recTime(call_k) = batInfo.(treat)(file_k).recTime;
                            cData.batNum{call_k} = batInfo.(treat)(file_k).batNum;
                            cData.daysOld(call_k) = round(days(batInfo.(treat)(file_k).recTime - cData.birthDates(strcmp(cData.batNums,cData.batNum{call_k}))));
                            cData.treatmentType(call_k) = batInfo.(treat)(file_k).treatmentType;
                            cData.file_call_pos(call_k,:) = batInfo.(treat)(file_k).callPos(file_call_k,:);
                            call_k = call_k + 1;
                        end
                    end
                end
                
            elseif strcmp(cData.expType,'piezo_recording')
                
                cData.loadWF = false;
                dataDir = 'C:\Users\phyllo\Documents\Maimon\acoustic_recording\';
                audio_base_dir = 'Z:\users\Maimon\acoustic_recording\audio\';
                recordingLogs = readtable([dataDir 'recording_logs.csv']);
                recordingLogs = recordingLogs(logical(recordingLogs.usable),:);
                treatmentGroups = readtable([dataDir 'bat_info.csv']);
                bat_ID_table_idx = contains(recordingLogs.Properties.VariableNames,'Bat_');
                nExp = size(recordingLogs,1);
                
                cData.batNums = cellfun(@num2str,num2cell(treatmentGroups.BatNum),'un',0);
                
                cData.nBats = length(cData.batNums);
                cData.birthDates = cell(1,cData.nBats);
                
                for b = 1:cData.nBats
                    dob = treatmentGroups.DOB(b);
                    cData.birthDates{b} = dob;
                end
                
                [cData.callWF, cData.batNum, cData.treatment, cData.fName, cData.fName_cut] = deal(cell(cData.maxCalls,1));
                [cData.loggerNum, cData.callLength, cData.daysOld] = deal(zeros(cData.maxCalls,1));
                cData.callPos = zeros(cData.maxCalls,2);
                cData.file_call_pos = zeros(cData.maxCalls,2);
                cData.expDay = datetime([],[],[]);
                
                date_str_format = 'mmddyyyy';
                cutFName = 'cut_call_data.mat';
                analysis_dir_name = 'Analyzed_auto';
                
                call_k = 1;
                
                for d = 1:nExp % iterate across all recording days
                    expDate = recordingLogs.Date(d);
                    exp_day_bat_nums = recordingLogs{d,bat_ID_table_idx};
                    exp_day_logger_nums = recordingLogs{d,find(bat_ID_table_idx)+1};
                    logger_bat_ID_idx = ~isnan(exp_day_bat_nums) & ~isnan(exp_day_logger_nums) & ~ismember(exp_day_logger_nums,str2double(recordingLogs.malfunction_loggers{1}));
                    
                    exp_day_logger_nums = exp_day_logger_nums(logger_bat_ID_idx);
                    exp_day_bat_nums = exp_day_bat_nums(logger_bat_ID_idx);
                    
                    dateStr = datestr(expDate,date_str_format);
                    audioDir = fullfile(audio_base_dir,dateStr,'audio\ch1\');
                    
                    cut_call_data = load(fullfile(audioDir,cutFName));
                    cut_call_data = cut_call_data.cut_call_data;
                    all_cut_call_files = dir(fullfile(audioDir,analysis_dir_name,'*Call*.mat'));
                    
                    assert(length(all_cut_call_files) == length(cut_call_data));
                    
                    AL_info = load(fullfile(audioDir,'AL_class_info.mat'));
                    
                    cut_call_data = cut_call_data(AL_info.usableIdx);
                    all_cut_call_files = all_cut_call_files(AL_info.usableIdx);
                    
                    predicted_loggerID = predict_AL_identity(AL_info);
                    noise_idx = cellfun(@length,predicted_loggerID) == 0;
                    
                    cut_call_data = cut_call_data(~noise_idx);
                    predicted_loggerID = predicted_loggerID(~noise_idx);
                    all_cut_call_files = all_cut_call_files(~noise_idx);
                    
                    identifiable_call_idx = cellfun(@length,predicted_loggerID) == 1;
                    usable_loggerID = [predicted_loggerID{identifiable_call_idx}];
                    loggerID = nan(1,sum(~noise_idx));
                    loggerID(identifiable_call_idx) = AL_info.logger_nums(usable_loggerID);
                    
                    if isempty(cut_call_data)
                        continue
                    end
                    
                    call_pos_expDay = vertcat(cut_call_data.corrected_callpos)/1e3; % convert to seconds
                    file_call_pos_expDay = vertcat(cut_call_data.callpos);
                    call_length_expDay = arrayfun(@(x) length(x.cut)/cData.fs,cut_call_data);
                    
                    for c = 1:length(cut_call_data)
                        cData.callLength(call_k) = call_length_expDay(c);
                        cData.callPos(call_k,:) = call_pos_expDay(c,:);
                        cData.file_call_pos(call_k,:) = file_call_pos_expDay(c,:);
                        cData.expDay(call_k,1) = expDate;
                        cData.batNum{call_k} = num2str(exp_day_bat_nums(exp_day_logger_nums == loggerID(c)));
                        cData.loggerNum(call_k) = loggerID(c);
                        cData.fName{call_k} = cut_call_data(c).fName;
                        cData.fName_cut{call_k} = fullfile(all_cut_call_files(c).folder,all_cut_call_files(c).name);
                        cData.callID(call_k) = cut_call_data(c).uniqueID;
                        
                        if ~isempty(cData.batNum{call_k})
                            b = strcmp(cData.batNums,cData.batNum{call_k});
                            cData.daysOld(call_k) = days(expDate - cData.birthDates{b});
                            cData.treatment{call_k} = treatmentGroups.Treatment{b};
                        else
                            cData.batNum{call_k} = 'multiple_bats_detected';
                            cData.treatment{call_k} = 'multiple_bats_detected';
                        end
                        call_k = call_k + 1;
                    end
                end
                cData.nCalls = call_k-1;
                
                
            else
                ME = MException('CallData:inputError','Experiment Type not recognized');
                throw(ME);
            end
            
            % if we initialized different data structures to have a length
            % of 'maxCalls' go ahead and shorten those to remove empty
            % elements
            callProperties = properties(cData)';
            for prop = callProperties
                if all(size(cData.(prop{:})) == [cData.maxCalls,1])
                    cData.(prop{:}) = cData.(prop{:})(1:cData.nCalls);
                elseif size(cData.(prop{:}),1) == cData.maxCalls && size(cData.(prop{:}),1) > 1
                    cData.(prop{:}) = cData.(prop{:})(1:cData.nCalls,:);
                end
            end
            
            cData.callID = reshape(cData.callID,[cData.nCalls 1]);
            
        end
        function varargout = subsref(cData,S)
            if length(S) == 2
                switch S(1).type
                    case '()'
                        nSubs = length(S(1).subs);
                        if ~rem(nSubs,2)
                            callIdx = true(cData.nCalls,1);
                            for idx = 1:2:nSubs
                                switch S(1).subs{idx}
                                    case 'cellInfo'
                                        callIdx = callIdx & strcmp(cData.cellInfo,S(1).subs{idx+1});
                                    case 'daysOld'
                                        daysOldIdx = false(cData.nCalls,1); %make an index assuming all 0 initially
                                        for d = S(1).subs{idx+1} %for the input that you're searching for
                                            daysOldIdx = daysOldIdx | round(cData.daysOld)==d; %make true if the day of indexed call is listed
                                        end
                                        callIdx = callIdx & daysOldIdx;
                                    case 'expDay'
                                        if length(S(1).subs{idx+1}) == 1
                                            callIdx = callIdx & (cData.expDay == S(1).subs{idx+1});
                                        else
                                            callIdx = callIdx & (cData.expDay >= S(1).subs{idx+1}(1) & cData.expDay < S(1).subs{idx+1}(2));
                                        end
                                    case 'batNum'
                                        multi_bat_idx = cellfun(@iscell,cData.batNum);
                                        batNumIdx = false(cData.nCalls,1);
                                        for bN = S(1).subs{idx+1}
                                            batNumIdx(~multi_bat_idx) = batNumIdx(~multi_bat_idx) | strcmp(cData.batNum(~multi_bat_idx),bN);
                                            batNumIdx(multi_bat_idx) = batNumIdx(multi_bat_idx) | cellfun(@(bNum) any(contains(bNum,bN)),cData.batNum(multi_bat_idx));
                                        end
                                        callIdx = callIdx & batNumIdx;
                                    case 'treatment'
                                        treatmentIdx = false(cData.nCalls,1);
                                        for b = S(1).subs{idx+1}
                                            treatmentIdx = treatmentIdx | strcmp(cData.treatment,b);
                                        end
                                        callIdx = callIdx & treatmentIdx;
                                    case 'batName'
                                        batNameIdx = false(cData.nCalls,1);
                                        for bN = S(1).subs{idx+1}
                                            batNameIdx = batNameIdx | ~cellfun(@isempty, strfind(cData.batName,bN)); %cellfun(@(x) ~isempty(strfind(x,bN)),cData.batName);
                                        end
                                        callIdx = callIdx & batNameIdx;
                                        %callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.batName);
                                    case 'batchNum'
                                        batchNumIdx = false(cData.nCalls,1);
                                        for bN = S(1).subs{idx+1}
                                            batchNumIdx = batchNumIdx | cData.batchNum == bN;
                                        end
                                        callIdx = callIdx & batchNumIdx;
                                    case 'callID'
                                        IDIdx = false(cData.nCalls,1);
                                        for id = S(1).subs{idx+1}'
                                            IDIdx = IDIdx | cData.callID == id;
                                        end
                                        callIdx = callIdx & IDIdx;
                                    case 'callNum'
                                        callNumIdx = false(cData.nCalls,1);
                                        for rN = S(1).subs{idx+1}
                                            callNumIdx = callNumIdx | cData.callNum==rN;
                                        end
                                        callIdx = callIdx & callNumIdx;
                                        %callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.callNum);
                                    case 'callTrigger'
                                        callIdx = callIdx & (cData.callTrigger >= S(1).subs{idx+1}(1) & cData.callTrigger <= S(1).subs{idx+1}(2));
                                    case 'callType'
                                        callTypeIdx = false(cData.nCalls,1);
                                        for cT = S(1).subs{idx+1}
                                            callTypeIdx = callTypeIdx | cellfun(@(x) ~isempty(strfind(x,cT)),cData.callType);
                                        end
                                        callIdx = callIdx & callTypeIdx;
                                        %callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.callType);
                                    case 'micType'
                                        callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.micType);
                                    case 'recNum'
                                        recNumIdx = false(cData.nCalls,1);
                                        for rN = S(1).subs{idx+1}
                                            recNumIdx = recNumIdx | cData.recNum==rN;
                                        end
                                        callIdx = callIdx & recNumIdx;
                                        %callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.recNum);
                                    case 'sessionID'
                                        sessionIDIdx = false(cData.nCalls,1);
                                        for sI = S(1).subs{idx+1}
                                            sessionIDIdx = sessionIDIdx | cData.sessionID==sI;
                                        end
                                        callIdx = callIdx & sessionIDIdx;
                                        %                                        callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.sessionID);
                                    case 'callLength'
                                        callLengthIdx = false(cData.nCalls,1);
                                        for cL = S(1).subs{idx+1}
                                            callLengthIdx = callLengthIdx | cData.callLength==cL;
                                        end
                                        callIdx = callIdx & callLengthIdx;
                                    case 'sessionType'
                                        sessionTypeIdx = false(cData.nCalls,1);
                                        for sT = S(1).subs{idx+1}
                                            sessionTypeIdx = sessionTypeIdx | strcmp(sT,cData.sessionType);
                                        end
                                        callIdx = callIdx & sessionTypeIdx;
                                        %callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.sessionType);
                                    case 'xrun'
                                        callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.xrun);
                                    otherwise
                                        disp('indexing variable not recognized')
                                        return
                                end
                            end
                            if iscell(cData.(S(2).subs)(callIdx)) && ~any(cellfun(@ischar,cData.(S(2).subs)(callIdx)))
                                try
                                    varargout = {vertcat(cData.(S(2).subs)(callIdx))};
                                catch
                                    try
                                        varargout = {[cData.(S(2).subs){callIdx}]};
                                    catch err
                                        display(err)
                                        return
                                    end
                                end
                                
                            else
                                varargout = {cData.(S(2).subs)(callIdx,:)};
                            end
                        else
                            disp('Indexing in VocalData must come in pairs');
                            return
                        end
                    otherwise
                        switch S(2).type
                            case '{}'
                                try
                                    varargout = {vertcat(cData.(S(1).subs){S(2).subs{:}})};
                                catch
                                    try
                                        varargout = {[cData.(S(1).subs){S(2).subs{:}}]};
                                    catch err
                                        display(err)
                                        return
                                    end
                                end
                            otherwise
                                try
                                    varargout = {builtin('subsref',cData,S)};
                                catch err
                                    switch err.message
                                        case 'Too many output arguments.'
                                            builtin('subsref',cData,S);
                                        otherwise
                                            display(err)
                                            return
                                    end
                                end
                        end
                end
            else
                try
                    varargout = {builtin('subsref',cData,S)};
                catch err
                    switch err.message
                        case 'Too many output arguments.'
                            builtin('subsref',cData,S);
                        otherwise
                            display(err)
                            return
                    end
                end
            end
            
            
        end
        function [interCallInterval, interDayIdx] = getICI(cData)
            cPos = cData.callPos;
            interCallInterval = [0; cPos(2:end,1) - cPos(1:end-1,2)]; % ICI(k) = time between call(k) and previous call
            if nargout > 1
                interDayIdx = abs(diff([cData.expDay(1); cData.expDay]))>duration(1,0,0) | interCallInterval<0;
            end
        end
        function inter_bat_ici = calculate_inter_bat_ICI(cData)
            bat_k = 1;
            inter_bat_ici = cell(1, length(cData.batNums));
            for bNum = cData.batNums
                callIdx = find(strcmp(cData.batNum,bNum{1}));
                cPos = cData.callPos(callIdx,:);
                ICI = [Inf; cPos(2:end,1) - cPos(1:end-1,2)];
                ICI = ICI(ICI > 0 & ICI < 60*60*12);
                iciIdx = ICI > 1;
                used_call_IDs = cData.callID(callIdx(iciIdx));
                try
                    S = substruct('()',{'callID',used_call_IDs,'batNum',bNum},'.','expDay');
                    callTimes = cData.subsref(S) + hours(12) + seconds(cPos(iciIdx,1));
                    
                catch
                    S = substruct('()',{'callID',used_call_IDs},'.','callID');
                    y = cData.subsref(S);
                    uniqueY = unique(y);
                    idx = find(histcounts(y,uniqueY)~=1);
                    discardCalls = false(1,length(idx));
                    rep_id_nums = nan(1,length(idx));
                    for k = 1:length(idx)
                        rep_id_nums(k) = unique(y(y == uniqueY(idx(k))));
                        S = substruct('()',{'callID',rep_id_nums(k)},'.','batNum');
                        rep_bat_num = cData.subsref(S);
                        nRep = length(rep_bat_num);
                        discardCalls(k) = nRep > 1 && length(unique(rep_bat_num)) ~= nRep;
                    end
                    
                    used_call_IDs = setdiff(used_call_IDs,rep_id_nums(discardCalls));
                    S = substruct('()',{'callID',used_call_IDs,'batNum',bNum},'.','callPos');
                    cPos = cData.subsref(S);
                    S = substruct('()',{'callID',used_call_IDs,'batNum',bNum},'.','expDay');
                    callTimes = cData.subsref(S) + hours(12) + seconds(cPos(:,1));
                    
                end
                
                callIdx = find(~strcmp(cData.batNum,bNum{1}));
                callTimes_nonbat = cData.expDay(callIdx) + hours(12) + seconds(cData.callPos(callIdx,1));
                
                inter_bat_ici{bat_k} = zeros(1,length(callTimes));
                for k = 1:length(callTimes)
                    allICI = callTimes_nonbat - callTimes(k);
                    allICI = allICI(allICI > 0 & allICI < seconds(60*60*12));
                    if ~isempty(allICI)
                        inter_bat_ici{bat_k}(k) = seconds(min(allICI));
                    else
                        inter_bat_ici{bat_k}(k) = NaN;
                    end
                end
                bat_k = bat_k + 1;
            end
        end
        function bout_call_nums = get_call_bout_nums(cData,callNum,boutSeparation)
            call_k = find(cData.callID == callNum);
            if ~isempty(cData.exp_date_num)
                session_idx = find(cData.exp_date_num == cData.exp_date_num(call_k));
            else
                session_idx = find(cData.expDay == cData.expDay(call_k));
            end
            cPos = cData.callPos(session_idx,:);
            session_inter_call_interval = [Inf; cPos(2:end,1) - cPos(1:end-1,2)];
            session_callID = cData.callID(session_idx);
            session_call_k = call_k - session_idx(1);
            idx = [session_call_k - find(session_inter_call_interval(session_call_k:-1:1) > boutSeparation,1,'first')+1:session_call_k ,...
                session_call_k+1:session_call_k+find(session_inter_call_interval(session_call_k+1:end) > boutSeparation,1,'first')-1];
            
            bout_call_nums = unique(session_callID(idx),'stable');
            if isempty(bout_call_nums)
                bout_call_nums = callNum;
            end
        end
        function [boutDuration, boutICI] = get_bout_duration(cData)
            boutSeparation = 1;
            
            cLength = cData.callLength;
            
            [interCallInterval, interDayIdx] = getICI(cData);
            interBoutIdx = [find(interCallInterval>boutSeparation); length(interCallInterval)];
            
            boutDuration = zeros(1,length(interBoutIdx)-1);
            boutICI = cell(1,length(interBoutIdx)-1);
            
            for bout_k = 1:length(interBoutIdx)-1
                last_call_in_bout = interBoutIdx(bout_k+1)-1;
                calls_in_bout = interBoutIdx(bout_k)+1:last_call_in_bout;
                if length(calls_in_bout) > 2
                    if any(interDayIdx(calls_in_bout))
                        last_call_in_bout = calls_in_bout(find(interDayIdx(calls_in_bout),1,'first'))-1;
                    end
                    calls_in_bout = interBoutIdx(bout_k)+1:last_call_in_bout;
                    boutDuration(bout_k) = sum(cLength(calls_in_bout)) + sum(interCallInterval(calls_in_bout));
                    boutICI{bout_k} = interCallInterval(calls_in_bout);
                    if boutDuration(bout_k) < 0
                        keyboard
                    end
                else
                    boutDuration(bout_k) = NaN;
                end
                
            end
            boutICI = vertcat(boutICI{:});
            boutDuration = boutDuration(boutDuration~=0 & ~isnan(boutDuration));
        end             
        function all_call_info = get_all_bhv(cData)
            
            call_bhv_dirs = dir(fullfile(cData.baseDirs{1},'bhv_data','call_info*.mat'));
            n_call_bhv_files = length(call_bhv_dirs);
            all_call_info = cell(1,n_call_bhv_files);
            for k = 1:n_call_bhv_files
                s = load(fullfile(call_bhv_dirs(k).folder,call_bhv_dirs(k).name));
                all_call_info{k} = s.call_info;
            end
            all_call_info = vertcat(all_call_info{:});
            
        end
        function oldestAge = get.oldestAge(obj)
            oldestAge = max(obj.daysOld);
        end
        function youngestAge = get.youngestAge(obj)
            youngestAge = min(obj.daysOld);
        end
        function yinParams = get.yinParams(obj)
            
            % Assemble parameter structure for yin algorithm
            wSize = round(obj.integrationWindow*obj.fs);
            overlapSize = round(obj.overlap*obj.fs);
            yinParams = struct('thresh',obj.thresh,'sr',obj.fs,'wsize',wSize,...
                'hop',wSize - overlapSize,'range',[1 obj.fs],...
                'bufsize',obj.fs,'APthresh',2.5,...
                'maxf0',obj.maxF0,'minf0',obj.minF0);
            
        end
        function featureParams = get.featureParams(obj)
            % Assemble parameter structure for windowed feature
            % calculation
            featureParams = struct('windowLength',obj.windowLength,'fs',obj.fs,...
                'overlap',obj.overlap,'cepstralFlag',obj.cepstralFlag,...
                'spCorrFlag',obj.spCorrFlag,'tapers',obj.tapers,...
                'maxF0',obj.maxF0,'minF0',obj.minF0);
        end
        function callWF = getCallWF(cData,call_k)
            if cData.loadWF
                callWF = cData.callWF{call_k};
            else
                callWF = loadCallWF_onTheFly(cData,call_k);
            end
        end
        function plotSpectrogram(cData,axisHandle,callInfo,specParams)
            
            if nargin < 4
                winSize = 5e-3;
                stepSize = 4.9e-3;
                specParams = struct('spec_win_size',round(winSize*cData.fs),...
                    'spec_overlap_size',round(stepSize*cData.fs),'spec_nfft',2^14,...
                    'spec_ylims',[0 60],'spec_caxis_factor',0.45);
            end
            
            if length(callInfo) == 1
                if cData.loadWF
                    data = cData.callWF{callInfo};
                else
                    data = loadCallWF_onTheFly(cData,callInfo);
                end
            else
                data = callInfo;
            end
            
            addpath('C:\Users\phyllo\Documents\GitHub\SoundAnalysisBats\')
            if ~isfield(specParams,'spec_caxis_factor')
                specParams.spec_caxis_factor = 0.75;
            end
            
            if length(data) <= specParams.spec_win_size
                return
            end
            
            [~,f,t,ps] = spectrogram(data,gausswin(specParams.spec_win_size),specParams.spec_overlap_size,specParams.spec_nfft,cData.fs,'yaxis');
            if length(t) <= 1
                return
            end
            f = f*1e-3; % frequency in kHz
            t = t*1e3;
            t = [0 t(end-1)];
            ps = 10*log10(ps);
            ps(isinf(ps)) = NaN;
            imagesc(axisHandle,t,f,ps)
            cmap = spec_cmap();
            colormap(cmap);
            ylim(specParams.spec_ylims);
            caxis([min(ps(:))*specParams.spec_caxis_factor max(ps(:))]);
            set(axisHandle,'YDir','normal')
            
        end
        function cData = manual_classify_call_type(cData,varargin)
            randomizeFlag = false;
            overWriteFlag = true;
            orig_rec_plot_win = 5;
            specParams = struct('spec_win_size',512,'spec_overlap_size',500,'spec_nfft',1024,'fs',cData.fs,'spec_ylims',[0 60],'spec_caxis_factor',0.75);
            if isempty(cData.manual_call_class)
                cData.manual_call_class = cell(1,cData.nCalls);
            end
            
            if ~isempty(varargin)
                start_call_k = varargin{1};
            else
                start_call_k = 1;
            end
            
            last_orig_wav_data = '';
            call_idx = start_call_k:cData.nCalls;
            
            if ~overWriteFlag
                call_idx = setdiff(call_idx,find(~cellfun(@isempty,cData.manual_call_class)));
            end
            
            if randomizeFlag
                call_idx = call_idx(randperm(length(call_idx)));
            end
            
            for call_k = call_idx
                data = getCallWF(cData,call_k);
                playback_fs = min(cData.fs,200e3);
                sound(data,playback_fs);
                
                origRec_fName = cData.fName{call_k};
                if iscell(origRec_fName)
                    origRec_fName = origRec_fName{1};
                end
                load_orig_wav_data = ~strcmp(origRec_fName,last_orig_wav_data);
                
                if load_orig_wav_data
                    dataFull = audioread(origRec_fName);
                    last_orig_wav_data = origRec_fName;
                end
                
                subplot(2,1,1)
                cla
                plotSpectrogram(cData,gca,specParams,data)
                subplot(2,1,2);
                hold on
                if load_orig_wav_data
                    cla
                    plot((1:length(dataFull))/cData.fs,dataFull,'k');
                end
                plot((cData.file_call_pos(call_k,1)+(0:length(data)-1))/cData.fs,data);
                
                if length(dataFull)/cData.fs > orig_rec_plot_win
                    xlim([cData.file_call_pos(call_k,1)/cData.fs - orig_rec_plot_win/2 cData.file_call_pos(call_k,1)/cData.fs + orig_rec_plot_win/2])
                end
                
                display([datestr(cData.expDay(call_k)) ' Call ' num2str(call_k)]);
                
                repeat = 1;
                repeat_k = 1;
                while repeat
                    class = input('Call type?','s');
                    if isempty(class)
                        pause(0.1);
                        if repeat_k < 3
                            sound(data,playback_fs/(2*repeat_k));
                        else
                            startIdx = max(1,cData.file_call_pos(call_k,1) - (orig_rec_plot_win/2)*cData.fs);
                            endIdx = min(length(dataFull),cData.file_call_pos(call_k,1) + (orig_rec_plot_win/2)*cData.fs);
                            sound(dataFull(startIdx:endIdx),playback_fs);
                            repeat_k = 1;
                        end
                        repeat_k = repeat_k + 1;
                        
                    else
                        if strcmp(class,'l')
                            cData.manual_call_class{call_k} = 'loud';
                            repeat = 0;
                        elseif strcmp(class,'c')
                            cData.manual_call_class{call_k} = 'connector';
                            repeat = 0;
                        elseif strcmp(class,'n')
                            cData.manual_call_class{call_k} = 'noise';
                            repeat = 0;
                        elseif strcmp(class,'m')
                            cData.manual_call_class{call_k} = 'mating';
                            repeat = 0;
                        elseif strcmp(class,'t')
                            cData.manual_call_class{call_k} = 'trill';
                            repeat = 0;
                        elseif strcmp(class,'o')
                            cData.manual_call_class{call_k} = 'other';
                            repeat = 0;
                        elseif strcmp(class,'g')
                            cData.manual_call_class{call_k} = 'grumble';
                            repeat = 0;
                        elseif strcmp(class,'i')
                            cData.manual_call_class{call_k} = 'interesting';
                            repeat = 0;
                        elseif strcmp(class,'stop')
                            return
                        elseif strcmp(class,'pause')
                            keyboard
                        end
                    end
                end
            end
            
        end
    end
end

function all_cut_call_data = get_cut_call_data(cData)


switch cData.expType
    
    case 'adult_operant'
        switch cData.exp_session_type
            case 'communication'
                cut_call_fnames = dir(fullfile(cData.baseDirs{1},'call_data','*cut_call_data.mat'));
            case 'operant'
                cut_call_fnames = dir(fullfile(cData.baseDirs{1},'call_data','*cut_call_data_operant*.mat'));
        end
        
    case 'adult'
        cut_call_fnames = dir(fullfile(cData.baseDirs{1},'call_data',['*cut_' cData.callEcho '_data.mat']));
        
    case 'juvenile'
        cut_call_fnames = cell(1,cData.nBats);
        for b = 1:cData.nBats
            cut_call_fnames{b} = dir(fullfile(cData.baseDirs{b},['bat' cData.batNums{b}],'**',['*cut_' cData.callEcho '_data.mat']));
        end
        cut_call_fnames = vertcat(cut_call_fnames{:});
        
end

all_cut_call_data = cell(1,length(cut_call_fnames));

for d = 1:length(cut_call_fnames) % iterate across all recording days
    
    s = load(fullfile(cut_call_fnames(d).folder,cut_call_fnames(d).name));
    cut_call_data = s.cut_call_data;
    
    if ~isempty(cut_call_data)
        
        if strcmp(cData.callEcho,'call')
            if all(islogical([cut_call_data.noise]))
                cut_call_data = cut_call_data(~[cut_call_data.noise]);
            else
                disp('non-logical values found in noise field of cut_call_data');
                keyboard;
            end
        elseif strcmp(cData.callEcho,'echo')
            [cut_call_data.noise] = deal(false);
        end
        
        
        if strcmp(cData.callEcho,'echo') && strcmp(cData.expType,'juvenile')
            call_info_fName = dir(fullfile(cut_call_fnames(d).folder,'juv_call_info*echo.mat'));
            if ~isempty(call_info_fName)
                call_info_fName = fullfile(call_info_fName.folder,call_info_fName.name);
                if ~exist(call_info_fName,'file')
                    continue
                end
            else
                continue
            end
            
            assert(all([cut_call_data.uniqueID] == [call_info.callID]));
            
            echo_calls = arrayfun(@(x) strcmp(x.echoCall,'juvEcho'),call_info);
            cut_call_data = cut_call_data(echo_calls);
        end
        
        all_cut_call_data{d} = cut_call_data;
    end
    
end

end

function callWF = loadCallWF_onTheFly(cData,call_k)

if any(strcmp(cData.expType,{'adult_operant','adult','juvenile'}))
    if ischar(cData.fName{call_k})
        callWF = audioread(cData.fName{call_k},cData.file_call_pos(call_k,:));
    else
        callWF1 = audioread(cData.fName{call_k}{1},[cData.file_call_pos(call_k,1) Inf]);
        callWF2 = audioread(cData.fName{call_k}{2},[1 cData.file_call_pos(call_k,2)]);
        callWF = [callWF1; callWF2];
    end
    if ~isempty(cData.compensation)
        [N,D] = rat(cData.compensation.fs/cData.fs);
        callWF = resample(callWF,D,N);
        callWF = filter(cData.compensation.irc,1,callWF);
    end
    
    if strcmp(cData.exp_session_type,'operant')
       [b,a] = butter(4,cData.EW_mic_HP_freq/(cData.fs/2),'high');
       callWF = filtfilt(b,a,callWF);
    end
    
elseif strcmp(cData.expType,'pratData')
    callWF = audioread(cData.fName{call_k});
    callWF = callWF(cData.file_call_pos(call_k,1):cData.file_call_pos(call_k,2));
elseif strcmp(cData.expType,'pratData')
    if ~isempty(cData.fName{call_k})
        callWF = load([cData.baseDirs cData.fName{call_k}]);
        callWF = callWF.convData';
    end
elseif strcmp(cData.expType,'pratData')
    callWF = load(cData.fName{call_k});
    try
        callWF = callWF.cut;
    catch
        callWF = callWF.finalcut;
    end
    callWF = reshape(callWF,numel(callWF),1);
    
elseif strcmp(cData.expType,'pratData')
    callWF = load(cData.fName{call_k});
    callWF = callWF.cut;
    callWF = reshape(callWF,numel(callWF),1);
elseif strcmp(cData.expType,'pratData')
    callWF = load(cData.fName_cut{call_k});
    callWF = callWF.cut;
    callWF = reshape(callWF,numel(callWF),1);
    
else
    disp('No functionality to load callWF for the experiment type');
    keyboard;
end

end

function [F0,ap] = calculate_yin_F0(callWF,P,interpFlag)

% Adapted from Vidush Mukund vidush_mukund@berkeley.edu
% March 13, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function builds acts as a wrapper to call the YIN Algorithm package;
% based on the script written of Yosef Prat
%   P.minf0:    Hz - minimum expected F0 (default: 30 Hz)
%   P.maxf0:    Hz - maximum expected F0 (default: SR/(4*dsratio))
%   P.thresh:   threshold (default: 0.1)
%   P.relfag:   if ~0, thresh is relative to min of difference function (default: 1)
%   P.hop:      samples - interval between estimates (default: 32)
%   P.range:    samples - range of samples ([start stop]) to process
%   P.bufsize:  samples - size of computation buffer (default: 10000)
%	P.sr:		Hz - sampling rate (usually taken from file header)
%	P.wsize:	samples - integration window size (defaut: SR/minf0)
%	P.lpf:		Hz - intial low-pass filtering (default: SR/4)
%	P.shift		0: shift symmetric, 1: shift right, -1: shift left (default: 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% call the yin function from the yin package
nSamples = length(callWF);
R = yin(callWF,P);

% F0 = 2.^(R.f0(round(P.wsize/2/P.hop)+1:end) + log2(440)); % convert from octaves to Hz
% ap = R.ap0(round(P.wsize/2/P.hop)+1:end);
F0 = 2.^(R.f0 + log2(440));
ap = R.ap0;

if interpFlag
    
    ap(ap ~= real(ap)) = abs(ap(ap ~= real(ap)));
    ap = smooth(ap,5);
    
    F0(ap > P.APthresh) = NaN;
    
    % any frequency values that are below the accepted range of fundamentals or
    % above the maximum fundamental frequency are replced with NaN
    F0(F0<=P.minf0 | F0>=P.maxf0 | F0<= 110) = NaN;
    
    % drop all NaN values (i.e. f0 values outside of the accepted range of
    % possible frequency values)
    ap = ap(~isnan(F0));
    F0 = F0(~isnan(F0));
    
    if nargin < 5
        N = length(F0);
    end
    
    if length(F0) > 1
        % perform a 1-D interpolation of the fundamenetal frequencies
        t=(1:length(F0)).*P.hop/P.sr + 0.5/P.minf0;
        tt=linspace(t(1),nSamples/P.sr,2*N+1);
        tt = tt(2:2:end-1);
        F0 = interp1(t,F0,tt);
        ap = interp1(t,ap,tt);
    elseif isempty(F0)
        F0 = NaN;
        ap = NaN;
    end
    idx = ~isnan(F0);
    F0 = F0(idx);
    ap = ap(idx);
end

if isempty(F0)
    F0 = NaN;
    ap = NaN;
end

F0 = smoothdata(F0,'movmean',2);
ap = smoothdata(ap,'movmean',2);

end
function F0 = calculate_spCorr_F0(callWF,p)
maxLag = round((1/p.minF0)*p.fs);
r = xcorr(callWF, maxLag, 'coeff');
F0 = spPitchCorr(r, p.fs, p.maxF0, p.minF0);
end
function [f0] = spPitchCorr(r, fs, mxf, mnf)
% NAME
%   spPitchCorr - Pitch Estimation via Auto-correlation Method
% SYNOPSIS
%   [f0] = spPitchCorr(r, fs)
% DESCRIPTION
%   Estimate pitch frequencies via Cepstral method
% INPUTS
%   r        (vector) of size (maxlag*2+1)x1 which contains Corr coefficients.
%             Use spCorr.m
%   fs       (scalar) the sampling frequency of the original signal
% OUTPUTS
%   f0       (scalar) the estimated pitch
% AUTHOR
%   Naotoshi Seo, April 2008
% SEE ALSO
%   spCorr.m
% search for maximum  between 2ms (=500Hz) and 20ms (=50Hz)
%adjusted parameters are 4000 Hz and 500 Hz
ms2=floor(fs/mxf); % 2ms
ms20=floor(fs/mnf); % 20ms
% half is just mirror for real signal
r = r(floor(length(r)/2):end);
[maxi,idx]=max(r(ms2:ms20));
f0 = fs/(ms2+idx-1);
end
function [weinerEntropy, spectralEntropy, centroid, energyEntropy, RMS,...
    pitchGoodness, cepstralF0, spCorrF0Win] = getCallFeatures(callWF,P)

L = length(callWF);
L_frame = P.windowLength*P.fs;
L_step = L_frame - P.overlap*P.fs;
nFrame = floor((L-(L_frame-L_step))/L_step);

if nFrame > 0
    [weinerEntropy, spectralEntropy, centroid, energyEntropy, RMS,...
        pitchGoodness, cepstralF0, spCorrF0Win] = deal(zeros(1,nFrame));
    for fr = 1:nFrame
        frameIdx = ((fr-1)*L_step +1):(((fr-1)*L_step)+L_frame);
        frame = callWF(frameIdx);
        [weinerEntropy(fr), spectralEntropy(fr), centroid(fr),...
            energyEntropy(fr), RMS(fr), pitchGoodness(fr), cepstralF0(fr),...
            spCorrF0Win(fr)] = getFeatures(frame,P);
    end
else
    [weinerEntropy, spectralEntropy, centroid, energyEntropy, RMS,...
        pitchGoodness, cepstralF0, spCorrF0Win] = deal(NaN);
end


end
function [weinerEntropy, spectralEntropy, centroid, energyEntropy, RMS,...
    pitchGoodness, cepstralF0, spCorrF0] = getFeatures(frame,P)

fs = P.fs;
tapers = P.tapers;
minF0 = P.minF0;
maxF0 = P.maxF0;

[AFFT, F] = afft(frame,fs);

%wiener entropy
weinerEntropy = exp(sum(log(AFFT)) ./ length(AFFT)) ./ mean(AFFT);

%spectral entropy
spectralEntropy = sentropy(AFFT);

%center of spectral mass (centroid)
centroid = sum( F' .* AFFT ) ./ sum(AFFT);

%Energy entropy
energyEntropy = sentropy(frame);

%RMS
RMS = sqrt(mean(frame.^2));

% Cepstral F0 pitch goodness
if P.cepstralFlag
    [pitchGoodness, cepstralF0] = cepstrumAnalysis(frame,tapers,fs);
else
    pitchGoodness = NaN;
    cepstralF0 = NaN;
end

if P.spCorrFlag % if we want to use the spCorr algorithm
    spCorrParams = struct('fs',fs,'minF0',minF0,'maxF0',maxF0);
    spCorrF0 = calculate_spCorr_F0(frame,spCorrParams);
else
    spCorrF0 = NaN;
end


end
function [AFFT, F] = afft(sig,fs,nfft)

if size(sig,1) == 1
    sig = sig';
end

L = size(sig,1);

if nargin < 3 || isempty(nfft)
    nfft = 2^nextpow2(L);
end

F = fs/2*linspace(0,1,nfft/2+1);
AFFT = fft(sig,nfft)./L;
AFFT = 2*abs(AFFT(1:nfft/2+1,:));

end
function ent = sentropy(v)

if size(v,2) == 1
    v = v';
end

v = v + abs(min(min(v(:)),0));
v = v./(ones(size(v,1),1)*sum(v));
v(v==0) = 1;
ent = -sum(v.*log(v));

end
function [pitchGoodness, f0] = cepstrumAnalysis(frame,tapers,fs)

frame = frame.*tapers(:,1);
tmp=rceps(frame);
[pitchGoodness, idx] = max(tmp(25:floor(end/2)));
f0 = (1/idx)*fs;

end

