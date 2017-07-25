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
        lesionedBats % which bats are lesioned (logical)
        deafenedBats % which bats are deafened (logical)
        nBats % number of bats
        cepstralFlag = false % perform cepstral analysis?
        spCorrFlag = false % perform spCorr analysis?
        batName %name of bat pairs in training
        callNum %call number (subset of recNum)
        recNum %recording number for the training session
        callTrigger %dt for autTrain
        callType %'callOnly' or 'callTrigger'
        sessionType %experiment type in the autoTrain
        sessionID %numerical counter for each sessionTye
        micType %which box/microphone used
        xrun %dropped sample number
                
        fs % sampling rate
        expType % String indicating which experiment we are dealing with
        
        callWF % call waveform
        expDay % experimental date (in datetime format)
        recTime % time of recording (in datetime format)
        daysOld % number of days from birth for this call
        batNum % bat ID number correspoding to this call
        fName % file path to call file
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
        
        maxCalls = 80000
        minF0 = 100 % in Hz
        maxF0 = 20e3 % in Hz
        
        windowLength = 5e-3 % in ms
        overlap = 4e-3 % in ms
        thresh = 0.1
        
        loadWF = true % Load in all call waveforms or not
    end
    properties (Dependent)
        oldestAge
        youngestAge
    end
    methods
        function cData = callData(fs, expType) % Function to initialize, load, and perform calculations
            cData.fs = fs; 
            cData.expType = expType;
            
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
            
            tapers = dpss(cData.windowLength*cData.fs,1.5);
            
            if ~cData.loadWF
               cData.callLength = zeros(cData.nCalls,1); 
            end
            
            % iterate through all calls
            for call_k = 1:cData.nCalls
                
                if cData.loadWF % if we already loaded all the data into this object
                    callWF = cData.callWF{call_k};
                else % if we need to load this data on the fly
                    if ~isempty(cData.fName{call_k})
                    callWF = loadCallWF_onTheFly(cData,call_k);
                    cData.callLength(call_k) = length(callWF)/fs;
                    end
                end
                
                % Assemble parameter structure for yin algorithm
                nSamples = length(callWF);
                wSize = round(cData.windowLength*cData.fs);
                overlap = round(cData.overlap*cData.fs);
                pYin = struct('thresh',cData.thresh,'sr',cData.fs,'wsize',wSize,...
                    'hop',wSize - overlap,'range',[1 nSamples],...
                    'bufsize',nSamples+2,'APthresh',2.5,...
                    'maxf0',cData.maxF0,'minf0',cData.minF0);
                
                
                [f0, ap0] = calculate_yin_F0(callWF,pYin);
                cData.yinF0Win{call_k} = f0;
                cData.ap0Win{call_k} = ap0;
                
                % take the F0 for the window of the call with the "best"
                % (i.e. lowest) aperiodicity)
                [~, idx] = min(ap0);
                cData.yinF0(call_k) = f0(idx);
                cData.ap0(call_k) = ap0(idx);
                
                
                % Assemble parameter structure for windowed feature
                % calculation
                pFeatures = struct('windowLength',cData.windowLength,'fs',cData.fs,...
                    'overlap',cData.overlap,'cepstralFlag',cData.cepstralFlag,'tapers',tapers);
                
                % calculate those features
                [cData.weinerEntropyWin{call_k},cData.spectralEntropyWin{call_k},...
                    cData.centroidWin{call_k},cData.energyEntropyWin{call_k},...
                    cData.RMSWin{call_k},cData.pitchGoodnessWin{call_k},...
                    cData.cepstralF0Win{call_k}] = getCallFeatures(callWF,pFeatures);
                
                % store average values of those features
                cData.weinerEntropy(call_k) = mean(cData.weinerEntropyWin{call_k});
                cData.spectralEntropy(call_k) = mean(cData.spectralEntropyWin{call_k});
                cData.centroid(call_k) = mean(cData.centroidWin{call_k});
                cData.energyEntropy(call_k) = mean(cData.energyEntropyWin{call_k});
                cData.RMS(call_k) = mean(cData.RMSWin{call_k});
                
                if cData.spCorrFlag % if we want to use the spCorr algorithm
                    cData.spCorrF0(call_k) = calculate_spCorr_F0(cData,call_k);
                end

                if cData.cepstralFlag % if we want to use the cepstral algorithm
                    [pg, idx] = max(cData.pitchGoodnessWin{call_k});
                    cData.pitchGoodness(call_k) = pg;
                    cData.cepstralF0(call_k) = cData.cepstralF0Win{call_k}(idx);
                end
            end
            
        end
        function n = numArgumentsFromSubscript(obj,~,~)
            n = numel(obj);
        end
        function cData = load_call_data(cData)
            
            switch cData.expType
                case 'ephys' % Maimon's juvenile ephys experiments
                    cData.baseDirs = {'Z:\users\Maimon\ephys\','E:\ephys\juvenile_recording\','E:\ephys\juvenile_recording\'};
                    cData.batNums = {'71319','71284','71173'};
                    cData.dateFormat = 'yyyyMMdd';
                    cData.birthDates = {datetime(2016,4,23),datetime(2016,09,24),datetime(2016,09,21)};
                    cData.nBats = length(cData.batNums);
                    
                    [cData.callWF, cData.batNum, cData.fName] = deal(cell(cData.maxCalls,1));
                    [cData.daysOld, cData.callLength] = deal(zeros(cData.maxCalls,1));
                    cData.callPos = zeros(cData.maxCalls,2);
                    cData.file_call_pos = zeros(cData.maxCalls,2);
                    cData.expDay = datetime([],[],[]);
                    
                    call_k = 1;
                    for b = 1:cData.nBats % iterate across all experimental bats
                        nlgDirs = dir([cData.baseDirs{b} 'bat' cData.batNums{b} filesep '*neurologger*']);
                        for d = 1:length(nlgDirs) % iterate across all recording days
                            audioDir = [cData.baseDirs{b} 'bat' cData.batNums{b} filesep nlgDirs(d).name filesep 'audio\ch1\'];
                            % get call data and time in recording
                            % (corrected for clock drift)
                            s = load([audioDir 'cut_call_data.mat']);
                            cut_call_data = s.cut_call_data;
                            
                            if ~isempty([cut_call_data.f_num])
                                
                                cutCalls = {cut_call_data.cut};
                                call_pos_expDay = vertcat(cut_call_data.corrected_callpos)/1e3; % convert to seconds
                                file_call_pos_expDay = vertcat(cut_call_data.callpos);
                                call_length_expDay = cellfun(@length, cutCalls)/cData.fs;
                                
                                expDatetime = datetime(nlgDirs(d).name(end-7:end),'InputFormat',cData.dateFormat);
                                
                                for c = 1:length(cutCalls)
                                    if ~cut_call_data(c).noise
                                        cData.callWF{call_k} = cutCalls{c};
                                        cData.callLength(call_k) = call_length_expDay(c);
                                        cData.callPos(call_k,:) = call_pos_expDay(c,:);
                                        cData.file_call_pos(call_k,:) = file_call_pos_expDay(c,:);
                                        cData.expDay(call_k) = expDatetime;
                                        cData.batNum{call_k} = cData.batNums{b};
                                        cData.daysOld(call_k) = days(expDatetime - cData.birthDates{b});
                                        cData.fName{call_k} = cut_call_data(c).fName;
                                        call_k = call_k + 1;
                                    end
                                end
                            end
                        end
                    end
                   
                    
                case 'lesions'
                    cData.batNums = {'71309','71306','71303','71305','71308'};
                    cData.nBats = length(cData.batNums);
                    cData.baseDirs = 'Z:\users\Eva\ES lesion recordings\';
                    cData.dateFormat = 'yyyyMMdd';
                    cData.lesionedBats = logical([1 1 0 0 0]);
                    cData.birthDates = {datetime(2016,11,27),datetime(2016,12,11),datetime(2016,10,25),datetime(2016,11,16),datetime(2016,10,19)};
                    data_var_name = 'finalcut';
                    treatmentTypes = {'lesioned','unmanipulated'};
                    data_dir_strs = {'lesion','unmanipulated control'};
                    batGroups = [1 1 2 2 2];
                    preceding_date_str = 'Box1';
                    
                    [cData.callWF, cData.batNum, cData.fName, cData.treatment] = deal(cell(cData.maxCalls,1));
                    [cData.daysOld, cData.callLength] = deal(zeros(cData.maxCalls,1));
                    cData.expDay = datetime([],[],[]);
                    
                    call_k = 1;
                    for t = 1:length(treatmentTypes)
                        batIdx = batGroups == t;
                        avg_birth_date = mean([cData.birthDates{batIdx}]);
                        analyzed_audio_dir = [cData.baseDirs data_dir_strs{t} filesep];
                        callFiles = dir([analyzed_audio_dir preceding_date_str '*Call*.mat']);
                        if ~isempty(callFiles)
                            for c = 1:length(callFiles)
                                cData.treatment{call_k} = treatmentTypes{t};
                                idx = strfind(callFiles(c).name,preceding_date_str) + length(preceding_date_str) + 1;
                                exp_date_str = callFiles(c).name(idx:idx+length(cData.dateFormat)-1);
                                data = load([analyzed_audio_dir callFiles(c).name],data_var_name);
                                cutCall = data.(data_var_name);
                                cData.callWF{call_k} = cutCall';
                                cData.callLength(call_k) = length(cutCall)/cData.fs;
                                
                                cData.expDay(call_k) = datetime(exp_date_str,'inputFormat',cData.dateFormat);
                                cData.fName{call_k} = [analyzed_audio_dir callFiles(c).name];
                                if isdatetime(avg_birth_date)
                                    cData.daysOld(call_k) = days(cData.expDay{call_k} - avg_birth_date);
                                else
                                    cData.daysOld(call_k) = NaN;
                                end
                                call_k = call_k + 1;
                            end
                        end
                    end
                    
                case 'deafened'
                    
                    cData.baseDirs = 'E:\deafened_recordings\all_recordings\';
                    cData.birthDates = {datetime(2016,9,14),datetime(2016,9,27),datetime(2016,8,19),datetime(2016,09,21),NaN,NaN};
                    cData.deafenedBats = logical([1 1 0 0 0 0]);
                    cData.batNums = {'71315','71354','65696','71353','1','2'};
                    treatmentTypes = {'deaf','saline','adult'};
                    data_dir_strs = {'deafened','saline control','adult control'};
                    batGroups = [1 1 2 2 3 3];
                    
                    cData.dateFormat = 'yyyyMMdd';
                    data_var_name = 'finalcut';
                    preceding_date_str = 'Box1';
                    
                    
                    [cData.callWF, cData.fName, cData.treatment] = deal(cell(cData.maxCalls,1));
                    [cData.callLength, cData.daysOld] = deal(zeros(cData.maxCalls,1));
                    cData.expDay = datetime([],[],[]);
                    
                    call_k = 1;
                    for t = 1:length(treatmentTypes)
                        batIdx = batGroups == t;
                        avg_birth_date = mean([cData.birthDates{batIdx}]);
                        analyzed_audio_dir = [cData.baseDirs data_dir_strs{t} filesep];
                        callFiles = dir([analyzed_audio_dir preceding_date_str '*Call*.mat']);
                        if ~isempty(callFiles)
                            for c = 1:length(callFiles)
                                cData.treatment{call_k} = treatmentTypes{t};
                                idx = strfind(callFiles(c).name,preceding_date_str) + length(preceding_date_str) + 1;
                                exp_date_str = callFiles(c).name(idx:idx+length(cData.dateFormat)-1);
                                data = load([analyzed_audio_dir callFiles(c).name],data_var_name);
                                cutCall = data.(data_var_name);
                                cData.expDay(call_k) = datetime(exp_date_str,'inputFormat',cData.dateFormat);
                                if isdatetime(avg_birth_date)
                                    cData.daysOld(call_k) = days(cData.expDay{call_k} - avg_birth_date);
                                else
                                    cData.daysOld(call_k) = NaN;
                                end
                                cData.callWF{call_k} = cutCall';
                                cData.callLength(call_k) = length(cutCall)/cData.fs;
                                cData.fName{call_k} = [analyzed_audio_dir callFiles(c).name];
                                call_k = call_k + 1;
                            end
                        end
                    end
                    
                case 'autoTrain'
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
                    
                case 'pratData'
                    cData.loadWF = false;
                    cData.baseDirs = 'E:\Yossi_vocalization_data\';
                    [fileInfo, batMetadata, batGroupIDs, allFileIDs_struct] = get_fileInfo_and_treatment_fileIDs(cData.baseDirs);
                    allFileIDs = allFileIDs_struct.con(1:1000);
                    cData.birthDates = mean([batMetadata.DOB{ismember(batMetadata.ID,batGroupIDs.con)}]);
                    treatmentTypes = {'isolated','control'};
                    
                    cData.nCalls = length(allFileIDs);
                    [cData.fName, cData.treatment] = deal(cell(cData.nCalls,1));
                    [cData.callLength, cData.daysOld] = deal(zeros(cData.nCalls,1));
                    cData.recTime = datetime(zeros(0,3));
                    
                    for call_k = 1:length(allFileIDs)
                        cData.isolated(call_k) = ismember(allFileIDs(call_k),allFileIDs_struct.iso);
                        cData.fName{call_k} = [cData.baseDirs fileInfo.fName{allFileIDs(call_k)}];
                        cData.recTime(call_k) = fileInfo.recTime(allFileIDs(call_k));
                        cData.daysOld(call_k) = day(cData.recTime(call_k) - cData.birthDates);
                    end
                otherwise
                    ME = MException('CallData:inputError','Experiment Type not recognized');
                    throw(ME);                    
            end
            
            % if we initialized different data structures to have a length
            % of 'maxCalls' go ahead and shorten those to remove empty
            % elements
            cData.nCalls = call_k-1;
            callProperties = properties(cData)';
            for prop = callProperties
                if all(size(cData.(prop{:})) == [cData.maxCalls,1])
                    cData.(prop{:}) = cData.(prop{:})(1:cData.nCalls);
                elseif size(cData.(prop{:}),1) == cData.maxCalls && size(cData.(prop{:}),1) > 1
                    cData.(prop{:}) = cData.(prop{:})(1:cData.nCalls,:);
                end
            end
            
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
                                            daysOldIdx = daysOldIdx | cData.daysOld==d; %make true if the day of indexed call is listed
                                        end
                                        callIdx = callIdx & daysOldIdx;
                                    case 'expDay'
                                        callIdx = callIdx & (cData.expDay >= S(1).subs{idx+1}(1) & cData.expDay < S(1).subs{idx+1}(2));
                                    case 'batNum'
                                        batNumIdx = false(cData.nCalls,1);
                                        for b = S(1).subs{idx+1}
                                            batNumIdx = batNumIdx | strcmp(cData.batNum,b);
                                        end
                                        callIdx = callIdx & batNumIdx;
                                    case 'treatment'
                                        callIdx = callIdx & strcmp(cData.treatment,S(1).subs{idx+1});
                                    case 'batName'
                                        batNameIdx = false(cData.nCalls,1);
                                        for bN = S(1).subs{idx+1}
                                            batNameIdx = batNameIdx | ~cellfun(@isempty, strfind(cData.batName,bN)); %cellfun(@(x) ~isempty(strfind(x,bN)),cData.batName);
                                        end
                                        callIdx = callIdx & batNameIdx;
                                        %callIdx = callIdx & cellfun(@(x) ~isempty(strfind(x,S(1).subs{idx+1})),cData.batName);                                       
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
                                        display('indexing variable not recognized')
                                        return
                                end
                            end
                            if iscell(cData.(S(2).subs)(callIdx)) && ~any(cellfun(@ischar,cData.(S(2).subs)(callIdx)))
                                try
                                    varargout = {vertcat(cData.(S(2).subs){callIdx})};
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
                            display('Indexing in VocalData must come in pairs');
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
        function oldestAge = get.oldestAge(obj)
            oldestAge = max(obj.daysOld);
        end
        function youngestAge = get.youngestAge(obj)
            youngestAge = min(obj.daysOld);
        end
    end
end

function callWF = loadCallWF_onTheFly(cData,call_k)

switch cData.expType
    case 'pratData'
        callWF = audioread([cData.baseDirs cData.fName{call_k}]);
    case 'autoTrain'
       if ~isempty(cData.fName{call_k})
            callWF = load([cData.baseDirs cData.fName{call_k}]);
            callWF = callWF.convData';
       end
            otherwise
        display('No functionality to load callWF for the experiment type');
        keyboard;
end

end

function [F0,ap] = calculate_yin_F0(callWF,P)

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
if isempty(F0)
    F0 = NaN;
    ap = NaN;
end

end
function F0 = calculate_spCorr_F0(cData,call_k)
maxLag = round(cData.windowLength*cData.fs);
r = xcorr(cData.callWF{call_k}, maxLag, 'coeff');
F0 = spPitchCorr(r, cData.fs, cData.maxF0, cData.minF0);
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
function [weinerEntropy, spectralEntropy, centroid, energyEntropy, RMS, pitchGoodness, cepstralF0] = getCallFeatures(callWF,P)

L = length(callWF);
L_frame = P.windowLength*P.fs;
L_step = L_frame - P.overlap*P.fs;
nFrame = floor((L-(L_frame-L_step))/L_step);

if nFrame > 0
    [weinerEntropy, spectralEntropy, centroid, energyEntropy, RMS, pitchGoodness, cepstralF0] = deal(zeros(1,nFrame));
    for fr = 1:nFrame
        frameIdx = ((fr-1)*L_step +1):(((fr-1)*L_step)+L_frame);
        frame = callWF(frameIdx);
        [weinerEntropy(fr), spectralEntropy(fr), centroid(fr), energyEntropy(fr), RMS(fr), pitchGoodness(fr), cepstralF0(fr)] = getFeatures(frame,P);
    end
else 
    [weinerEntropy, spectralEntropy, centroid, energyEntropy, RMS, pitchGoodness, cepstralF0] = deal(NaN);
end


end
function [weinerEntropy, spectralEntropy, centroid, energyEntropy, RMS, pitchGoodness, cepstralF0] = getFeatures(frame,P)

fs = P.fs;
tapers = P.tapers;

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
