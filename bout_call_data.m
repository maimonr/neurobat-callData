classdef bout_call_data < ephysData
    properties
        inter_bout_thresh = 1 % seconds
        min_calls_in_bout = 2
        min_bout_length = 0.4 % seconds
        fs = 250e3
        envelope_window_size = 4e-3 % seconds
        lp_filt_order = 4;
        lp_filt_stopband_freq = 300;
        downsample_factor = 1e3
        frequency_resolution = 1 % Hz
        boutSize = 25 % seconds
        spikeOffset = [-1 1];
        feature_threshold = 1e-3
        rhythm_freq_band = [3 8];
        rhythm_slidin_win_s = 0.5
        rhythm_sliding_win_step_s = 0.04
        rhythm_lag_segment_bounds = linspace(3,10,4);
        rhythm_coef_thresh = 0.05
        rhythm_min_seq_length_s = 0.25
        boot_reps = 1e3
        rSpec_type = 'mt'
        freq_res = 5
        frequencyLims = [1 50]
        callType = 'call'
        
        rSpec
        nBouts
        boutCalls
        boutEnv
        batNum
        cellNums
        cellInfo
        callPos
        file_call_pos
        daysOld
        expDay
        boutSpikes
        sortingQuality
        rhythmSeqs
        bout_spike_rhythm_phase
        treatment
        raw_data_filter
    end
    
    properties (Dependent)
        feature_fs
        freq
        nfft
        time_half_bandwidth
        spikeRange
        rhythm_corr_lags
        boutLength
    end
    
    methods
        
        function bc = bout_call_data(cData,vdCall)
            
            bc@ephysData(cData.expType);
            
            debugFlag = false;
            bc.batNums = cData.batNums;
            bc.fs = cData.fs;
            
            if strcmp(bc.expType,'adult_operant') && strcmp(bc.callType,'operant')
                [b,a] = butter(4,200/(bc.fs/2),'high');
                bc.raw_data_filter = struct('b',b,'a',a);
            end
            
            [b,a] = butter(bc.lp_filt_order,bc.lp_filt_stopband_freq/(bc.fs/2),'low');
            feature_filter_coefs = [b;a];
            
            [b_bp,a_bp] = butter(bc.lp_filt_order,[bc.rhythm_freq_band(1)/(bc.feature_fs/2),bc.rhythm_freq_band(2)/(bc.feature_fs/2)],'bandpass');
            [b_lp,a_lp] = butter(bc.lp_filt_order,bc.rhythm_freq_band(2)/(bc.feature_fs/2),'low');
            rhythm_filter_coefs = struct('lp',[b_lp;a_lp],'bp',[b_bp;a_bp]);
            
            nBat = length(bc.batNums);
            
            [rhythmSpec,boutCalls,features,rhythmSeqs,expDays] = deal(cell(1,nBat));
            
            for b = 1:nBat
                batNum = bc.batNums(b);
                if any(isnan(cData('batNum',batNum).daysOld) | cData('batNum',batNum).daysOld==0)
                    expDays{b} = cData('batNum',batNum).expDay;
                    daysOffset = (24*60*60)*days(expDays{b} - expDays{b}(1));
                else
                    daysOffset = (24*60*60)*cData('batNum',batNum).daysOld;
                end
                callPos = cData('batNum',batNum).callPos + daysOffset;
                callIDs = cData('batNum',batNum).callID;
                file_callPos = cData('batNum',batNum).file_call_pos;
                interCallInterval = [Inf; callPos(2:end,1) - callPos(1:end-1,2)]; % ICI(k) = time between call(k) and previous call
                interCallInterval(interCallInterval<0) = abs(interCallInterval(interCallInterval<0));
                
                interBoutIdx = [find(interCallInterval>bc.inter_bout_thresh); length(interCallInterval)];
                nCalls = length(interBoutIdx)-1;
                rhythmSpec{b} = nan(nCalls,(bc.nfft/2)+1);
                boutCalls{b} = cell(1,nCalls);
                features{b} = cell(nCalls,1);
                rhythmSeqs{b} = cell(nCalls,1);
                k = 1;
                lastProgress = 0;
                for bout_k = 1:nCalls
                    
                    calls_in_bout = callIDs(interBoutIdx(bout_k):interBoutIdx(bout_k+1)-1);
                    call_fNames = cData('callID',calls_in_bout).fName;
                    n_calls_in_bout = length(calls_in_bout);
                    
                    bout_batNums = cData('callID',calls_in_bout).batNum;
                    
                    if n_calls_in_bout<bc.min_calls_in_bout || any(cellfun(@iscell,bout_batNums)) || length(unique(bout_batNums))>1
                        success = false;
                    else
                        bout_file_callpos = [file_callPos(interBoutIdx(bout_k),1) file_callPos(interBoutIdx(bout_k+1)-1,2)];
                        [success,boutCallWF] = get_bout_WF(bc,call_fNames,bout_file_callpos);
                    end
                    
                    if success
                        if debugFlag
                            sound(boutCallWF,bc.fs)
                        end
                        boutCalls{b}{bout_k} = calls_in_bout;
                        features{b}{bout_k} = get_bout_feature(bc,boutCallWF,feature_filter_coefs);
                        rhythmSpec{b}(bout_k,:) = get_bout_rSpec(bc,features{b}{bout_k});
                        rhythmSeqs{b}{bout_k} = get_rhythm_bout_timestamps(bc,features{b}{bout_k},rhythm_filter_coefs);
                    end
                    [k,lastProgress] = updateProgress(nCalls,k,lastProgress);
                end
                rhythmSpec{b} = rhythmSpec{b}(~cellfun(@isempty,boutCalls{b}),:);
                features{b} = features{b}(~cellfun(@isempty,boutCalls{b}));
                rhythmSeqs{b} = rhythmSeqs{b}(~cellfun(@isempty,boutCalls{b}));
                expDays{b} = expDays{b}(~cellfun(@isempty,boutCalls{b}));
                boutCalls{b} = boutCalls{b}(~cellfun(@isempty,boutCalls{b}));
            end
            bc.expDay = vertcat(expDays{:});
            bc.rSpec = vertcat(rhythmSpec{:});
            bc.boutCalls = horzcat(boutCalls{:})';
            bc.boutEnv = vertcat(features{:});
            bc.rhythmSeqs = vertcat(rhythmSeqs{:});
            bc.nBouts = size(bc.rSpec,1);
            
            bNum = cell(1,length(bc.batNums));
            for b = 1:length(bc.batNums)
                bNum{b} = repmat(bc.batNums(b),1,length(boutCalls{b}));
            end
            
            bc.batNum = horzcat(bNum{:})';
            
            bc.daysOld = cellfun(@(x) cData('callID',x(1)).daysOld,bc.boutCalls,'un',0);
            bc.daysOld = cellfun(@(x) x(1),bc.daysOld);
            
            bout_call_pos = cellfun(@(x) [cData('callID',x(1)).callPos cData('callID',x(end)).callPos],bc.boutCalls,'un',0);
            bc.callPos = cellfun(@(x) [x(1) x(end)]',bout_call_pos,'un',0);
            
            bout_file_call_pos = cellfun(@(x) [cData('callID',x(1)).file_call_pos cData('callID',x(end)).file_call_pos],bc.boutCalls,'un',0);
            bc.file_call_pos = cellfun(@(x) [x(1) x(end)]',bout_file_call_pos,'un',0);
            
            if nargin > 1
                bc.cellNums = cellfun(@(x) find(cellfun(@(y) ismember(x(1),y),vdCall.callNum)),bc.boutCalls,'un',0);
                bc.cellInfo = cellfun(@(x) vdCall.cellInfo(cellfun(@(y) ismember(x(1),y),vdCall.callNum)),bc.boutCalls,'un',0);
                bc.sortingQuality = cellfun(@(x) vdCall.sortingQuality(cellfun(@(y) ismember(x(1),y),vdCall.callNum)),bc.boutCalls,'un',0);
                
                
                bc.boutSpikes = cell(bc.nBouts,1);
                bc.bout_spike_rhythm_phase = cell(bc.nBouts,1);
                for bout_k = 1:bc.nBouts
                    bc.boutSpikes{bout_k} = get_bout_spike_times(bc,vdCall,bout_k);
                    bc.bout_spike_rhythm_phase{bout_k} = get_spike_rhythm_phase(bc,bout_k,rhythm_filter_coefs);
                end
            end
        end
        
        function n = numArgumentsFromSubscript(~,~,~)
            n = 1;
        end
        
        function varargout = subsref(bc,S)
            if length(S) == 2
                switch S(1).type
                    case '()'
                        nSubs = length(S(1).subs);
                        if ~rem(nSubs,2)
                            boutIdx = true(bc.nBouts,1);
                            for idx = 1:2:nSubs
                                switch S(1).subs{idx}
                                    case 'cellInfo'
                                        boutIdx = boutIdx & cellfun(@(x) any(strcmp(x,S(1).subs{idx+1})),bc.cellInfo);
                                    case 'daysOld'
                                        daysOldIdx = false(1,bc.nBouts);
                                        for d = S(1).subs{idx+1}
                                            daysOldIdx = daysOldIdx | bc.daysOld==d;
                                        end
                                        boutIdx = boutIdx & daysOldIdx;
                                    case 'expDay'
                                        if length(S(1).subs{idx+1}) == 1
                                            boutIdx = boutIdx & (bc.expDay == S(1).subs{idx+1});
                                        else
                                            boutIdx = boutIdx & (bc.expDay >= S(1).subs{idx+1}(1) & bc.expDay < S(1).subs{idx+1}(2));
                                        end
                                    case 'batNum'
                                        batNumIdx = false(bc.nBouts,1);
                                        for b = S(1).subs{idx+1}
                                            batNumIdx = batNumIdx | strcmp(bc.batNum,b);
                                        end
                                        boutIdx = boutIdx & batNumIdx;
                                    otherwise
                                        display('indexing variable not recognized')
                                        return
                                end
                            end
                            if iscell(bc.(S(2).subs)(boutIdx)) && ~any(cellfun(@ischar,bc.(S(2).subs)(boutIdx)))
                                try
                                    varargout = {bc.(S(2).subs)(boutIdx)};
                                    if any(strcmp(S(1).subs,'cellInfo'))
                                        idx = find(strcmp(S(1).subs,'cellInfo'));
                                        bout_cell_info = bc.cellInfo(boutIdx);
                                        cellIdx = cellfun(@(cellInfo) find(strcmp(cellInfo,S(1).subs{idx+1})),bout_cell_info,'un',0);
                                        varargout = {cellfun(@(x,idx) x(idx),varargout{1},cellIdx)};
                                    end
                                catch err
                                    display(err)
                                    return
                                end
                            else
                                if size(bc.(S(2).subs),1) == length(boutIdx)
                                    varargout = {bc.(S(2).subs)(boutIdx,:)};
                                elseif size(bc.(S(2).subs),2) == length(boutIdx)
                                    varargout = {bc.(S(2).subs)(:,boutIdx)};
                                else
                                    disp('Indexing mismatch');
                                    return
                                end
                            end
                        else
                            disp('Indexing in VocalData must come in pairs');
                            return
                        end
                    otherwise
                        switch S(2).type
                            case '{}'
                                try
                                    varargout = {vertcat(bc.(S(1).subs){S(2).subs{:}})};
                                catch
                                    try
                                        varargout = {[bc.(S(1).subs){S(2).subs{:}}]};
                                    catch err
                                        display(err)
                                        return
                                    end
                                end
                            otherwise
                                try
                                    varargout = {builtin('subsref',bc,S)};
                                catch err
                                    switch err.message
                                        case 'Too many output arguments.'
                                            builtin('subsref',bc,S);
                                        otherwise
                                            display(err)
                                            return
                                    end
                                end
                        end
                end
            else
                try
                    varargout = {builtin('subsref',bc,S)};
                catch err
                    switch err.message
                        case 'Too many output arguments.'
                            builtin('subsref',bc,S);
                        otherwise
                            display(err)
                            return
                    end
                end
            end
            
            
        end
        
        function [p,spike_bout_rhythm_phase] = spike_phase_shuffle_test(bc,cell_k)
            addpath('C:\Users\phyllo\Documents\MATLAB\circstat\')
            [b_bp,a_bp] = butter(bc.lp_filt_order,[bc.rhythm_freq_band(1)/(bc.feature_fs/2),bc.rhythm_freq_band(2)/(bc.feature_fs/2)],'bandpass');
            [b_lp,a_lp] = butter(bc.lp_filt_order,bc.rhythm_freq_band(2)/(bc.feature_fs/2),'low');
            filterCoefs = struct('lp',[b_lp;a_lp],'bp',[b_bp;a_bp]);
            
            bout_ks = find(cellfun(@(x) ismember(cell_k,x),bc.cellNums))';
            bout_ks = bout_ks(~cellfun(@isempty,bc.rhythmSeqs(bout_ks)));
            if ~isempty(bout_ks)
                spike_bout_rhythm_phase_shuffled = cell(1,length(bout_ks));
                spike_bout_rhythm_phase = cellfun(@(x,y) x(y==cell_k),bc.bout_spike_rhythm_phase(bout_ks),bc.cellNums(bout_ks));
                spike_bout_rhythm_phase = [spike_bout_rhythm_phase{:}];
                
                for k = 1:length(bout_ks)
                    cell_num = find(bc.cellNums{bout_ks(k)} == cell_k);
                    spike_bout_rhythm_phase_shuffled{k} = get_spike_rhythm_phase_shuffled(bc,bout_ks(k),cell_num,filterCoefs);
                end
                spike_bout_rhythm_phase_shuffled = [spike_bout_rhythm_phase_shuffled{:}];
                if ~isempty(spike_bout_rhythm_phase_shuffled)
                    shuffle_resultant_length = circ_r(spike_bout_rhythm_phase_shuffled');
                    true_resultant_length = circ_r(spike_bout_rhythm_phase');
                    p = sum(shuffle_resultant_length >= true_resultant_length)/bc.boot_reps;
                    
                else
                    p = NaN;
                end
                
            else
                p = NaN;
                spike_bout_rhythm_phase = NaN;
            end
            
            
        end
        
        function boutCallWF = get_bout_WF(bc,cData,bout_k,offset)
            
            if nargin < 4
                offset = 0;
            end
            
            calls_in_bout = bc.boutCalls{bout_k};
            
            call_fNames = cData('callID',calls_in_bout).fName;
            
            sample_call_offset = bc.fs*offset*[-1; 1];
            [~,boutCallWF] = get_bout_WF(bc,call_fNames,bc.file_call_pos{bout_k}+sample_call_offset);
            
            if ~isempty(cData.compensation)
                [N,D] = rat(cData.compensation.fs/cData.fs);
                boutCallWF = resample(boutCallWF,D,N);
                boutCallWF = filter(cData.compensation.irc,1,boutCallWF);
            end
            
            
        end
        
        function [success,boutCallWF,sample_offset] = get_bout_WF_manual(bc,cData,call_fNames,file_callPos,offset)
            sample_call_offset = bc.fs*offset;
            sample_call_offset = reshape(sample_call_offset,size(file_callPos));
            file_callPos_offset = file_callPos + sample_call_offset;
            
            [success,boutCallWF,sample_offset] = get_bout_WF(bc,call_fNames,file_callPos_offset);
            
            if success && ~isempty(cData.compensation)
                [N,D] = rat(cData.compensation.fs/cData.fs);
                boutCallWF = resample(boutCallWF,D,N);
                boutCallWF = filter(cData.compensation.irc,1,boutCallWF);
            end
        end
        
        function rhythmic_seq_timestamps = get_rhythm_ts_manual(bc,bout_k)
            [b_bp,a_bp] = butter(bc.lp_filt_order,[bc.rhythm_freq_band(1)/(bc.feature_fs/2),bc.rhythm_freq_band(2)/(bc.feature_fs/2)],'bandpass');
            [b_lp,a_lp] = butter(bc.lp_filt_order,bc.rhythm_freq_band(2)/(bc.feature_fs/2),'low');
            rhythm_filter_coefs = struct('lp',[b_lp;a_lp],'bp',[b_bp;a_bp]);
            
            rhythmic_seq_timestamps = get_rhythm_bout_timestamps(bc,bc.boutEnv{bout_k},rhythm_filter_coefs);
        end
        
        function boutSpikes = get_bout_spike_times(bc,vdCall,bout_k,cell_ks)
            
            if nargin < 4
                cell_ks = bc.cellNums{bout_k};
            end
            
            call_spike_range = bc.callPos{bout_k}(1) + bc.spikeRange;
            nCells = length(cell_ks);
            boutSpikes = cell(1,nCells);
            
            for k = 1:nCells
                timestamps = getSpikes(vdCall,cell_ks(k));
                timestamps = 1e-3*timestamps;
                boutSpikes{k} = inRange(timestamps,call_spike_range) - bc.callPos{bout_k}(1);
            end
            
        end
        
        function [callWF, callPos_cat, nSample] = concatenate_bouts(bc,cData,call_ks,silenceRange)
            
            nCall = length(call_ks);
            callWF = [];
            callPos_cat = [];
            nSample = zeros(1,nCall+1);
            k = 1;
            for call_k = call_ks
                
                bout_call_pos = vertcat(cData('callID',bc.boutCalls{call_k}).file_call_pos);
                bout_call_pos = bout_call_pos - bout_call_pos(1) + 1 + nSample(k);
                if any(bout_call_pos < 0)
                    continue
                end
                
                boutCallWF = bc.get_bout_WF(cData,call_k);
                n_noise_sample = randi(silenceRange*bc.fs);
                
                callWF = [callWF boutCallWF zeros(1,n_noise_sample)];
                callPos_cat = [callPos_cat; bout_call_pos];
                k = k + 1;
                nSample(k) = length(callWF);
            end
        end
        
        function feature_fs = get.feature_fs(obj)
            feature_fs = round(obj.fs/obj.downsample_factor);
        end
        
        function nfft = get.nfft(obj)
            switch obj.rSpec_type
                case 'mt'
                    nfft = 2^(1+nextpow2(obj.inter_bout_thresh*obj.feature_fs));
                case 'avg_spectrogram'
                    nfft = 2^(1+nextpow2(obj.inter_bout_thresh*obj.feature_fs));
                case 'persistence'
                    nfft =  2*(diff(obj.frequencyLims)+1)*(obj.freq_res-1);
            end
        end
        
        function freq = get.freq(obj)
            switch obj.rSpec_type
                case 'mt'
                    freq = linspace(0,obj.feature_fs/2,(obj.nfft/2)+1);
                case 'persistence'
                    freq = linspace(obj.frequencyLims(1),obj.frequencyLims(2),(obj.nfft/2)+1);
            end
            
        end
        
        function time_half_bandwidth = get.time_half_bandwidth(obj)
            time_half_bandwidth = obj.frequency_resolution*(1/obj.feature_fs)*(obj.boutSize*obj.feature_fs);
        end
        
        function spikeRange = get.spikeRange(obj)
            spikeRange = [obj.spikeOffset(1) obj.boutSize+obj.spikeOffset(2)];
        end
        
        function boutLength = get.boutLength(obj)
            
            boutLength =  diff([obj.callPos{:}]);
            
        end
        
    end
end

function rSpec = get_bout_rSpec(bc,feature)

debugFlag = false;
feature_orig = feature;
if size(feature,1) ~= 1
    feature = feature';
end

try
    thresh = bc.feature_threshold;
    feature(feature>thresh) = 1;
    feature(feature<thresh) = 0;
    
    feature = feature - mean(feature);
    
    switch bc.rSpec_type
        case 'mt'
            padSize = ceil(((bc.boutSize*bc.feature_fs)-length(feature))/2);
            padFeature = [zeros(1,padSize) feature zeros(1,padSize)];
            
            if rem(length(feature),2) ~= 0
                padFeature = padFeature(1:end-1);
            end
            
            if length(padFeature)/bc.feature_fs ~= bc.boutSize
                keyboard
            end
            
            if rem(length(feature),2) ~= 0
                padFeature = padFeature(1:end-1);
            end
            rSpec = pmtm(padFeature,bc.time_half_bandwidth,bc.nfft,bc.feature_fs);
        case 'persistence'
            
            [p,f_persistence,pwr] = pspectrum(feature,bc.feature_fs,'persistence','FrequencyLimits',bc.frequencyLims,'FrequencyResolution',5,'Leakage',0.75);
            [X,Y] = meshgrid(f_persistence,pwr);
            rSpec = mean(X.*Y.*p*(length(feature)/bc.feature_fs),1);
            
        case 'avg_spectrogram'
            
            if length(feature)/bc.feature_fs < bc.inter_bout_thresh
                winSize = length(feature);
                rSpec = periodogram(feature,gausswin(winSize),bc.nfft,bc.feature_fs);
            else
                winSize = round(bc.inter_bout_thresh*bc.feature_fs);
                [~,~,~,ps] = spectrogram(feature,gausswin(winSize),round(0.9*winSize),bc.nfft,bc.feature_fs);
                rSpec = mean(ps,2);
            end
            
    end
catch err
    disp(err)
    keyboard
end

if debugFlag
    clf(figure(1))
    clf(figure(2))
    max_f = 50;
    min_f = 1;
    [p,f_persistence,pwr] = pspectrum(feature,bc.feature_fs,'persistence','FrequencyLimits',[min_f max_f],'FrequencyResolution',5,'Leakage',0.75);
    [X,Y] = meshgrid(f_persistence,pwr);
    
    figure(1);
    subplot(3,1,1)
    plot((1:length(feature))/bc.feature_fs,feature)
    
    subplot(3,1,2)
    plot((1:length(feature_orig))/bc.feature_fs,feature_orig)
    
    subplot(3,1,3)
    hold on
    plot(bc.freq(bc.freq<max_f & bc.freq>min_f),((rSpec(bc.freq<max_f & bc.freq>min_f))))
    %
    if length(feature)/bc.feature_fs < bc.inter_bout_thresh
        winSize = length(feature);
        [pxx,f] = periodogram(feature,gausswin(winSize),bc.nfft,bc.feature_fs);
    else
        winSize = round(bc.inter_bout_thresh*bc.feature_fs);
        [~,f,~,ps] = spectrogram(feature,gausswin(winSize),round(0.9*winSize),bc.nfft,bc.feature_fs);
        pxx = mean(ps,2);
    end
    
    
    plot(f(f<max_f & f>min_f),((pxx(f<max_f & f>min_f))))
    xlim([0 max_f])
    
    
    figure(2);
    pspectrum(feature,bc.feature_fs,'persistence','FrequencyLimits',[min_f max_f],'FrequencyResolution',5,'Leakage',0.75);
    
    debug_pause = input('pause?');
    
    if debug_pause
        keyboard;
    end
    
    clf(figure(1))
    clf(figure(2))
    
end

end

function feature = get_bout_feature(bc,boutCallWF,filterCoefs)

senv = envelope(boutCallWF,bc.envelope_window_size*bc.fs,'rms');
senv = filtfilt(filterCoefs(1,:),filterCoefs(2,:),senv);
feature = downsample(senv,bc.downsample_factor);

end

function [success,boutCallWF,sample_offset] = get_bout_WF(bc,call_fNames,file_callPos)

nCalls = length(call_fNames);
if any(cellfun(@iscell,call_fNames))
    multi_file_fnames = call_fNames{cellfun(@iscell,call_fNames)};
    call_fNames{cellfun(@iscell,call_fNames)} = multi_file_fnames{2};
end

fNames_in_bout = unique(call_fNames);

file_nums = cellfun(@(fname) strsplit(fname,filesep),fNames_in_bout,'un',0);
file_nums = cellfun(@(fname) strsplit(fname{end}(1:strfind(fname{end},'.')-1),'_'),file_nums,'un',0);
file_nums = cellfun(@(fname) str2double(fname{end}),file_nums);
[~,idx] = sort(file_nums,'ascend');
[audio_dir,~,ext] = fileparts(fNames_in_bout{idx(1)});

all_fnames = dir(fullfile(audio_dir,['*' ext]));
all_fnums = arrayfun(@(fname) strsplit(fname.name,'_'),all_fnames,'un',0);
all_fnums = cellfun(@(fname) strsplit(fname{end}(1:strfind(fname{end},'.')-1),'_'),all_fnums,'un',0);
all_fnums = cellfun(@(fname) str2double(fname{end}),all_fnums);

if file_callPos(1) < 1
    first_bout_fnum = file_nums(idx(1));
    
    preceding_file_idx = all_fnums == first_bout_fnum-1;
    preceding_fname = fullfile(audio_dir,all_fnames(preceding_file_idx).name);
    info = audioinfo(preceding_fname);
    fNames_in_bout = [preceding_fname; fNames_in_bout];
    file_callPos(1) = file_callPos(1) + info.TotalSamples;
end

last_bout_fname = fNames_in_bout{idx(end)};
info = audioinfo(last_bout_fname);
if file_callPos(2) > info.TotalSamples
    last_bout_fnum = file_nums(idx(end));
    
    next_file_idx = all_fnums == last_bout_fnum+1;
    next_fname = fullfile(audio_dir,all_fnames(next_file_idx).name);
    fNames_in_bout = [fNames_in_bout; next_fname];
    file_callPos(2) = file_callPos(2) - info.TotalSamples;
end

n_wav_files = length(fNames_in_bout);
file_callPos_by_file = zeros(n_wav_files,2);

boutCallWF = [];
nSamp_per_file = zeros(1,n_wav_files+1);
for f = 1:n_wav_files
    info = audioinfo(fNames_in_bout{f});
    nSamp_per_file(f+1) = info.TotalSamples;
    if f == 1
        file_callPos_by_file(f,:) = file_callPos;
    else
        file_callPos_by_file(f-1,2) = nSamp_per_file(f);
        file_callPos_by_file(f,:) = [1 file_callPos(2)];
    end
end
nSamp_per_file = cumsum(nSamp_per_file);
boutCallPos = [file_callPos(1) nSamp_per_file(n_wav_files)+file_callPos(2)];
if abs(diff(boutCallPos))/bc.fs < bc.min_bout_length || abs(diff(boutCallPos))/bc.fs > bc.boutSize
    success = false;
    boutCallWF = [];
    sample_offset = [];
    return
end

for f = 1:n_wav_files
    callData = audioread(fNames_in_bout{f},file_callPos_by_file(f,:));
    callData = reshape(callData,1,numel(callData));
    boutCallWF = [boutCallWF callData];
end

if ~isempty(bc.raw_data_filter)
    boutCallWF = filtfilt(bc.raw_data_filter.b,bc.raw_data_filter.a,boutCallWF);
end

sample_offset = zeros(1,nCalls);

for k = 1:nCalls
    sample_offset(k) = nSamp_per_file(strcmp(fNames_in_bout,call_fNames{k}));
end

success = true;

end

function rhythmic_seq_timestamps = get_rhythm_bout_timestamps(bc,feature,filter_coefs)

debugFlag = false;
slinding_win_size = round(bc.feature_fs*bc.rhythm_slidin_win_s);
sliding_win_step_size = round(bc.feature_fs*bc.rhythm_sliding_win_step_s);
lags = (-(slinding_win_size-1):(slinding_win_size-1))/bc.feature_fs;

freqs = 1./lags;

lagIdx = false(length(bc.rhythm_lag_segment_bounds)-1,length(lags));

for k = 1:length(bc.rhythm_lag_segment_bounds)-1
    lagIdx(k,:) = freqs>bc.rhythm_lag_segment_bounds(k) & freqs<bc.rhythm_lag_segment_bounds(k+1);
end
%%
padFeature = [zeros(1,slinding_win_size) feature zeros(1,slinding_win_size)];
paddedEnv = filtfilt(filter_coefs.lp(1,:),filter_coefs.lp(2,:),padFeature');
% paddedEnv = zscore(paddedEnv);
paddedEnv = diff(paddedEnv);
% thresh = bc.feature_threshold*std(diff(zscore(feature)));
thresh = bc.feature_threshold;
paddedEnv(paddedEnv>thresh) = 1;
paddedEnv(paddedEnv<thresh) = 0;

idx = slidingWin(length(paddedEnv),slinding_win_size,slinding_win_size-sliding_win_step_size);
r = zeros(size(idx,1),(2*size(idx,2))-1);
for k = 1:size(idx,1)
    r(k,:) = xcorr(paddedEnv(idx(k,:)).*kaiser(size(idx,2)),'unbiased');
end

pad_t = ((0:length(paddedEnv)-1)/bc.feature_fs)-bc.rhythm_slidin_win_s;
sliding_win_t = mean(pad_t(idx),2);
bout_idx = find(sliding_win_t>0 & sliding_win_t<length(feature)/bc.feature_fs);
bout_idx = [bout_idx(1)-1; bout_idx; bout_idx(end)+1];

sliding_win_t = sliding_win_t(bout_idx);
r = abs(r(bout_idx,:));

lag_segment_coefs = zeros(size(lagIdx,1),size(r,1));
for k = 1:size(lagIdx,1)
    lag_segment_coefs(k,:) = abs(movmean(mean((r(:,lagIdx(k,:))),2),round(bc.rhythm_slidin_win_s/bc.rhythm_sliding_win_step_s)));
end

seqs_indices = lag_segment_coefs > bc.rhythm_coef_thresh;
seqs_indices = remove_short_segments(seqs_indices,round(bc.rhythm_min_seq_length_s/bc.rhythm_sliding_win_step_s));
seqs_indices = any(seqs_indices);

if any(seqs_indices)
    rhythmic_seq_timestamps = [sliding_win_t(diff([0 seqs_indices])==1) sliding_win_t(diff([seqs_indices 0])==-1)];
    rhythmic_seq_timestamps(rhythmic_seq_timestamps<0) = 0;
else
    rhythmic_seq_timestamps = [];
end

if debugFlag
    
    t = (0:length(feature)-1)/bc.feature_fs;
    subplot(3,1,1);
    cla
    hold on
    plot([0 max(t)],repmat(bc.rhythm_coef_thresh,1,2),'--r')
    plot(sliding_win_t,lag_segment_coefs)
    plot(sliding_win_t(seqs_indices),0.5*ones(1,sum(seqs_indices)),'kx')
    
    xlim([0 max(t)])
    %     ylim([0 0.5])
    
    subplot(3,1,2);
    cla
    hold on
    plot(t,feature)
    plot(t,filtfilt(filter_coefs.lp(1,:),filter_coefs.lp(2,:),feature))
    xlim([0 max(t)])
    
    h = hilbert(filtfilt(filter_coefs.bp(1,:),filter_coefs.bp(2,:),feature));
    instant_phase = angle(h);
    subplot(3,1,3);
    plot(t,instant_phase)
    xlim([0 max(t)])
    
    debug_pause = input('?');
    
    if debug_pause
        keyboard
    end
end
%%
end

function spike_bout_rhythm_phase = get_spike_rhythm_phase(bc,bout_k,filter_coefs)
debugFlag = false;

boutSpikes = bc.boutSpikes{bout_k};
spike_bout_rhythm_phase = cell(1,length(boutSpikes));

rhythmic_seq_timestamps = bc.rhythmSeqs{bout_k};

if ~isempty(rhythmic_seq_timestamps)
    t = (0:length(bc.boutEnv{bout_k})-1)/bc.feature_fs;
    instant_phase = get_bout_phase(bc.boutEnv{bout_k},filter_coefs,t);
    for cell_k = 1:length(boutSpikes)
        spike_idx = false(1,length(boutSpikes{cell_k}));
        for k = 1:size(rhythmic_seq_timestamps,1)
            [~,idx] = inRange(boutSpikes{cell_k},rhythmic_seq_timestamps(k,:));
            spike_idx = spike_idx | idx;
        end
        spikes = boutSpikes{cell_k}(spike_idx);
        spike_bout_rhythm_phase{cell_k} = nan(1,sum(spike_idx));
        for spike_k = 1:sum(spike_idx)
            [~,spike_t_idx] = min(abs(t - spikes(spike_k)));
            spike_bout_rhythm_phase{cell_k}(spike_k) = instant_phase(spike_t_idx);
        end
    end
    if debugFlag
        keyboard;
    end
end
end

function spike_bout_rhythm_phase = get_spike_rhythm_phase_shuffled(bc,bout_k,cell_k,filter_coefs)

boutSpikes = bc.boutSpikes{bout_k}{cell_k};

timestamps = bc.spikeRange(1):1e-3:bc.spikeRange(2);

rhythmic_seq_timestamps = bc.rhythmSeqs{bout_k};

if ~isempty(rhythmic_seq_timestamps)
    t = (0:length(bc.boutEnv{bout_k})-1)/bc.feature_fs;
    instant_phase = get_bout_phase(bc.boutEnv{bout_k},filter_coefs,t);
    spike_idx = false(1,length(boutSpikes));
    timestamps_idx = false(1,length(timestamps));
    for k = 1:size(rhythmic_seq_timestamps,1)
        [~,idx] = inRange(boutSpikes,rhythmic_seq_timestamps(k,:));
        spike_idx = spike_idx | idx;
        
        [~,idx] = inRange(timestamps,rhythmic_seq_timestamps(k,:));
        timestamps_idx = timestamps_idx | idx;
    end
    boutSpikes = boutSpikes(spike_idx);
    timestamps = timestamps(timestamps_idx);
    
    nSpike = length(boutSpikes);
    spike_bout_rhythm_phase = zeros(bc.boot_reps,nSpike);
    for b = 1:bc.boot_reps
        spikes = timestamps(randperm(length(timestamps),nSpike));
        for spike_k = 1:sum(spike_idx)
            [~,spike_t_idx] = min(abs(t - spikes(spike_k)));
            spike_bout_rhythm_phase(b,spike_k) = instant_phase(spike_t_idx);
        end
        
    end
else
    spike_bout_rhythm_phase = [];
end
end

function instant_phase = get_bout_phase(feature,filter_coefs,t)

filt_env = filtfilt(filter_coefs.bp(1,:),filter_coefs.bp(2,:),feature);
%     h = hilbert(filt_env);
%     instant_phase = angle(h);
f = fit(t',filt_env','sin1');
instant_phase = mod(((f.b1*t)+f.c1),2*pi)-pi;

end

function seqs = remove_short_segments(x,length_thresh_samples)

seqs = x;
dx = diff([zeros(size(x,1),1) x zeros(size(x,1),1)],[],2);
for k = 1:size(x,1)
    for d = find(dx(k,:)==1)
        seq_samples = find(dx(k,d:end)==-1,1,'first')-1;
        if isempty(seq_samples)
            seq_samples = length(d:size(dx,2));
        end
        
        if seq_samples < length_thresh_samples
            end_seq_idx = d - 1 + seq_samples;
            seqs(k,d:end_seq_idx) = false;
        end
    end
end

end

function [k,lastProgress] = updateProgress(nCalls,k,lastProgress)

progress = 100*(k/nCalls);

if mod(progress,10) < mod(lastProgress,10)
    fprintf('%d %% of calls processed\n',round(progress));
end

lastProgress = progress;

k = k + 1;

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