% Coded by J.L-C.
% Newencode Analytics
% www.newencode.com
% April-September, 2022
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% CLEAR THE WORKSPACE
clear ; clc
close all
dt1 = datetime;
% ##########################################################################
% ##########################  EDIT STARTS  #################################
% ##########################################################################
% ReadOptionParams = 0 read params from the header specs of this script
% ReadOptionParams = 1 load params from file (mat file)
ReadOptionParams = 0 ;

TaskName = 'BEES';

% folder where you have the raw data (mff files)
path_source   = '/Volumes/WORKBENCH/LISADATA/BEES_ERR_mff files';

% folder where you have the channel location file
chanlocfolder = '/Volumes/WORKBENCH/LISADATA/BCDLAB_FTAG_beta002/acquisition_layout_information';
chanlocfile   = 'GSN-HydroCel-129.sfp'; % channel location file

% parameter file location
paramFilePath = '/Volumes/WORKBENCH/LISADATA/BCDLAB_FTAG_beta002/code';
paramFileName = 'LISA_BEBES_1.mat';

%
% Subjects' spreadsheet (be sure you are using the right spreadsheet per session (pre/post) )
%
spsheet = '/Volumes/WORKBENCH/LISADATA/BEES_Infant_erptrain_Counterbalance.xlsx';

typeFields = {'code'};
onsetTags = { 'ct02' 'it01' } ;

removelinenoise = false; % true or false; remove line noise (60 Hz in the US)
linefreq = 60;

HP_type    = 1;   % 1: HIGH-PASS FILTER; 2: POLYNOMIAL DETRENDING
hpf_cutoff = 0.05; % high-pass filter cutoff in Hz (e.g. 0.1)
lpf_cutoff = 30;  % low-pass filter cutoff in Hz (e.g. 30)

removeouterband = false;
chouterband ={ 'E17' 'E43' 'E48' 'E49' 'E56' 'E63' 'E68'...
        'E73' 'E81' 'E88' 'E94' 'E99' 'E107' 'E113'...
        'E119' 'E120' 'E123' 'E124' 'E125' 'E126'...
        'E127' 'E128' } ;
% max percentage (e.g. 10, 20, 50) of bad channels allowed (to be interpolated)
maxperbadch2int = 50 ;

rawnchan        = 129; % number of channels of raw data
workingchannels = 129; % correct number of (final) channels
workingsrate    = 500; % correct sample rate

%   8 = 'E8'
%   9 = 'E9'
%  14 = 'E14'
%  15 = 'E15'
%  20 = 'E21'
%  21 = 'E22'
%  24 = 'E25'
InterpolbadChans = true;
eblnkchan = [8 9 14 15 18 20 21 24]; % eye movement channels

%
% if you want this script to set your Matlab path -automatically- by adding the required
% packages then make loadPackagesIntoMatlab = true and set the location of
% the BCDLAB_FTAG folder
%
loadPackagesIntoMatlab = true;
BCDLAB_FTAG_fullpath = '/Volumes/WORKBENCH/LISADATA/BCDLAB_FTAG_beta002';

% ##########################################################################
% ##########################   EDIT ENDS  ##################################
% ##########################################################################

if loadPackagesIntoMatlab

        % SET FOLDERS FOR HAPPE AND EEGLAB PATHS-----------------------------------
        fprintf('Preparing HAPPE...\n') ;

        % SET HAPPE AND EEGLAB PATHS USING THE RUNNING SCRIPT----------------------
        happeDir  = fileparts(which(mfilename('fullpath'))) ;
        happeDir  = strrep(happeDir, [filesep 'Code' filesep 'Newencode'],'');
        eeglabDir = [happeDir filesep 'Packages' filesep 'eeglab2021.0'] ;

        % ADD HAPPE AND REQUIRED FOLDERS TO PATH-----------------------------------
        addpath([happeDir filesep 'acquisition_layout_information'], ...
                [happeDir filesep 'Code'], ...
                [happeDir filesep 'Code' filesep 'Newencode'], ...
                [happeDir filesep 'scripts'], ...
                [happeDir filesep 'scripts' filesep 'UI_scripts'], ...
                [happeDir filesep 'scripts' filesep 'pipeline_scripts'], ...
                eeglabDir, genpath([eeglabDir filesep 'functions'])) ;
        rmpath(genpath([eeglabDir filesep 'functions' filesep 'octavefunc'])) ;

        % ADD EEGLAB FOLDERS TO PATH-----------------------------------------------
        pluginDir = dir([eeglabDir filesep 'plugins']) ;
        pluginDir = strcat(eeglabDir, filesep, 'plugins', filesep, {pluginDir.name}, ';') ;
        addpath([pluginDir{:}]) ;

        % ADD CLEANLINE FOLDERS TO PATH--------------------------------------------
        if exist('cleanline', 'file')
                cleanlineDir = which('eegplugin_cleanline.m') ;
                cleanlineDir = cleanlineDir(1:strfind(cleanlineDir, 'eegplugin_cleanline.m')-1) ;
                addpath(genpath(cleanlineDir)) ;
        else; error('Please make sure cleanline is on your path') ;
        end

        % ADD ERPLAB FOLDERS TO PATH-----------------------------------------------
        if exist('erplab_default_values', 'file')
                erplabDir = which('eegplugin_erplab.m') ;
                erplabDir = erplabDir(1:strfind(erplabDir, 'eegplugin_erplab.m')-1) ;
                addpath(genpath(erplabDir)) ;
        else; error('Please make sure ERPLAB is on your path') ;
        end
end

switch ReadOptionParams
        case 0
                params.lowDensity         = 0;
                params.paradigm.task      = 1 ;
                params.paradigm.ERP.on    = 0;
                params.paradigm.onsetTags = onsetTags ;
                params.QCfreqs            = [ 2:3:8 12 20 30 45 70];
                params.loadInfo.inputFormat   = 3 ;
                params.loadInfo.layout        = [ 2 128] ;
                params.loadInfo.correctEGI    = 1 ;
                params.loadInfo.chanlocs.inc  = 1;
                params.loadInfo.chanlocs.file = fullfile(chanlocfolder, chanlocfile);
                params.loadInfo.typeFields    = typeFields ;
                params.chans.subset        =  []  ;
                params.chans.IDs           = chouterband;
                params.lineNoise.freq      = linefreq ;
                params.lineNoise.neighbors = [] ;
                params.lineNoise.legacy    = 0 ;
                params.downsample          = workingsrate;
                params.filt.butter         = [] ;
                params.badChans.rej        = InterpolbadChans ;
                params.badChans.legacy     = [] ;
                params.wavelet.legacy      = [] ;
                params.segment.on          = [];
                params.segment.interp      = InterpolbadChans ;
                params.segRej.on           = [] ;
                params.reref.on            = 1 ;
                params.reref.chan          = 'Cz'  ;
                params.reref.average       = 1 ;
                params.reref.subset        = { 'Cz' } ;
                params.outputFormat        = [];
                params.vis.enabled         = 0 ;
                params.vis.min             = 0 ;
                params.vis.max             = 0 ;
                params.vis.toPlot          = [ ] ;
                params.HAPPEver            = '2_2_0' ;
        case 1
                % LOAD THE PARAMETER FILE--------------------------------------
                disp('Loading parameters...') ;
                load(fullfile(paramFilePath,paramFileName)) ;
                disp('Parameters loaded.')
        otherwise

end

% SHOW PARAMETERS----------------------------------------------------------
listParams(params) ;
nc_showParamsStruct(params) ; %Javier's
lnMeans  = [] ;
wavMeans = [] ;

% ------------------------------------------------------------------------
%
% load subjects' spreadsheet
%
Tsubj = readtable(spsheet);
Tsubj = Tsubj(:, 1:5); % only first 4 columns are valid
Tsubj = standardizeMissing(Tsubj,{''});
Tsubj = rmmissing(Tsubj);
Tsubj = renamevars(Tsubj,1:5, ...
        ["ID","AGE","CB","Gender","Observation"]);
Tsubj = convertvars(Tsubj,["ID","AGE","CB","Gender","Observation"],'categorical');

% PREPARE DATASET LIST TO BE PROCESSED-------------------------------------
a = dir(path_source);
b = regexpi({a.name}','.*(?=\.)','match');
b = [b{:}];
% only get subjects from IDs at the spreadsheet
b = b(contains(b, append("_",string([Tsubj.ID])',"_"),'IgnoreCase',true));
% only get subjects from the right task
b = b(contains(b,TaskName,'IgnoreCase',true));

if ~isempty(b)
        fileArrayAux = strtrim(unique(b)); % potential filenames to be loaded
        fileArray = cell(1,10);
        hh = 1;
        for qq = 1:length(fileArrayAux) % check the existency of the file signal1.bin
                if exist(fullfile(path_source,[fileArrayAux{qq} '.mff'],'signal1.bin'), 'file') == 2
                        fileArray{hh} = fileArrayAux{qq};
                        hh = hh + 1;
                else
                        %file does not exist.
                end
        end
        fileArray = fileArray(~cellfun('isempty',fileArray)); % final array of filenames
        nfile = length(fileArray);
else
        nfile = 0;
end

% INITIALIZE QUALITY REPORT METRICS----------------------------------------
fprintf('Initializing report metrics...\n') ;
% DATA QUALITY METRICS: create a variable holding all the names of each
% metric and a variable to hold the metrics for each file.
dataQCnames = {'File_Length_in_Seconds', 'Number_User-Selected_Chans', ...
        'Number_Good_Chans_Selected', 'Percent_Good_Chans_Selected', 'Bad_Chan_IDs', ...
        'Percent_Var_Retained_Post-Wav', 'Chans_Interpolated',...
        'Filename','Age_from_file','True_age','Obs_from_file','True_Obs',...
        'CB_from_file','True_CB','ct02_count','it01_count','Amount_Samples',...
        'Splitting_sample','QC'};

dataQC     = cell(length(fileArray), length(dataQCnames)) ;

% #########################################################################
%                               DATASET LOOP
% #########################################################################
propbadch = maxperbadch2int/100; % bad chan
for kdata=1:nfile
%         try
                % init variables for output table
                Age_from_file = NaN;
                true_AGE      = NaN;
                Obs_from_file = 'null';
                true_Obs      = 'null';
                true_CB       = NaN;
                ct02_count = NaN;
                it01_count = NaN;
                totalpnts = NaN;
                point2split = NaN;

                % current filename (raw data .mff)
                sname = fileArray{kdata};

                % CB from filename
                get_CB = contains(sname,'cb2','IgnoreCase',true);
                if get_CB
                        CB_from_file = 2;
                else
                        CB_from_file = 1;
                end

                % Examples of file names
                % BEES_121_12_1_erptrain_20180427_112505.mff

                % subject ID from filename
                subj_info = regexpi(sname,'_(\d{3}+)_(\d{1}+)*[_]*(\d{1}+)*','tokens'); % get ID and age (if any)
                subj_info = subj_info{:};
                subj_ID   = subj_info{1}; % ID from filename

                % get age if any
                if length(subj_info)>=2
                        if isempty(subj_info{2})
                                subj_AGE = 'not_declared';
                        else
                                subj_AGE  = subj_info{2}; % age from filename
                        end
                else
                        subj_AGE = 'not_declared';
                end

                % get observation if any
                if length(subj_info)>=3
                        if isempty(subj_info{3})
                                subj_observa = 'not_declared';
                        else
                                subj_observa  = subj_info{3}; % age from filename
                        end
                else
                        subj_observa = 'not_declared';
                end

                % CB from spreadsheet
                cb_value = Tsubj.CB(Tsubj.ID==subj_ID);
                if size(cb_value,1)>1 && strcmpi(subj_observa,'not_declared')
                        error('dataset has missing information in its file name (e.g. "observation" value)')
                elseif size(cb_value,1)>1 && ~strcmpi(subj_observa,'not_declared')
                        cb_value = Tsubj.CB(Tsubj.ID==subj_ID & Tsubj.Observation==subj_observa);
                        if isempty(cb_value)
                                error('dataset''s filename does not match information in spreadsheet (e.g. "observation" value)')
                        end
                end
                true_CB  = str2num(char(cb_value));

                % AGE (calculation in months) from spreadsheet
                age_value = Tsubj.AGE(Tsubj.ID==subj_ID);

                if size(age_value,1)>1 && strcmpi(subj_observa,'not_declared')
                        error('dataset has missing information in its file name (e.g. "observation" value)')
                elseif size(age_value,1)>1 && ~strcmpi(subj_observa,'not_declared')
                        age_value = Tsubj.AGE(Tsubj.ID==subj_ID & Tsubj.Observation==subj_observa);
                        if isempty(age_value)
                                error('dataset''s filename does not match information in spreadsheet (e.g. "observation" value)')
                        end
                end
                true_AGE = regexpi(char(age_value), '(\d+)\s*Y[,]\s*(\d+)\s*M[,]\s*(\d+)\s*D','tokens');
                true_AGE = true_AGE{:};
                true_AGE = str2double(true_AGE);
                true_AGE = round(12*true_AGE(1) + true_AGE(2) + true_AGE(3)/30);

                if strcmpi(subj_AGE, 'not_declared')
                        subj_AGE = num2str(true_AGE);
                else
                        Age_from_file = str2num(subj_AGE);
                end

                % get Observation from spreadsheet
                obs_value = Tsubj.Observation(Tsubj.ID==subj_ID);
                if size(obs_value,1)>1 && strcmpi(subj_observa,'not_declared')
                        error('dataset has missing information in its file name (e.g. "observation" value)')
                elseif size(obs_value,1)>1 && ~strcmpi(subj_observa,'not_declared')
                        obs_value = Tsubj.Observation(Tsubj.ID==subj_ID & Tsubj.Observation==subj_observa);

                        if isempty(obs_value)
                                error('dataset''s filename does not match information in spreadsheet (e.g. "observation" value)')
                        end
                end

                % observation (visit the lab: 1,2,3,...)
                true_Obs  = char(obs_value);
                if strcmpi(subj_observa, 'not_declared')
                        subj_observa  = true_Obs;
                        Obs_from_file =  'not_declared';
                else
                        Obs_from_file = subj_observa;
                end

                %###############################################################
                % -----------------------LOAD RAW DATA-------------------------
                %###############################################################
                fprintf('\n%s\n\n', repmat('#',1, 150));
                fprintf('Loading %s.mff... Please wait\n', sname);
                fprintf('\n%s\n', repmat('#',1, 150));
                EEG = pop_mffimport({fullfile(path_source, sprintf('%s.mff', sname))},params.loadInfo.typeFields,0,0);

                % CHECK CHANNEL NUMBER (CHCK 1)--------------------------------------
                if EEG.nbchan~=rawnchan
                        error('wrong channel number')
                end

                % Report file length in seconds'
                dataQC{kdata, 1} = EEG.xmax ;

                % LOAD CHANNEL LOC INFO----------------------------------------------
                EEG = pop_chanedit(EEG, 'load', {params.loadInfo.chanlocs.file ...
                        'filetype' 'autodetect'}) ;
                % VALIDATE
                EEG = eeg_checkset(EEG) ;

                % #############################################################
                % REMOVE OUTER BAND OF ELECTRODES------------------------------------
                % #############################################################
                if removeouterband
                        chanIDs = setdiff({EEG.chanlocs.labels}, params.chans.IDs);
                        EEG     = pop_select(EEG, 'channel', chanIDs) ;
                        EEG.setname = [EEG.setname '_cs'] ;
                        EEG = eeg_checkset( EEG );
                end

                chanIDs = {EEG.chanlocs(1:end-1).labels}; % skip Cz

                % Report how many channels are present in the dataset as a
                % data metric.
                dataQC{kdata,2} = size(chanIDs,2);

                % CHECK CHANNEL NUMBER (CHCK 2)--------------------------------------
                if EEG.nbchan~=workingchannels
                        error('wrong channel number')
                end

                % #############################################################
                % #               DELETE IDLE EEG SEGMENT (ERPLAB)
                % #############################################################
                try
                        EEG = nc_eegtrim(EEG, 1000, 1000, onsetTags);
                catch
                        fprintf('\nOops!...There is not enough samples to keep the pre-stimulation window...\n');
                end

                % #############################################################
                % #         apply high-pass filter or polynomial detrending
                % #############################################################
                switch HP_type % skip last channel Cz
                        case 1
                                % HIGH-PASS FILTER hpf_cutoff Hz (ERPLAB)-----------------------------------
                                fprintf('Highpass filtering the data...\n');
                                EEG  = pop_basicfilter( EEG,  1:EEG.nbchan-1 , 'Boundary', 'boundary', 'Cutoff', hpf_cutoff,...
                                        'Design', 'butter', 'Filter', 'highpass',...
                                        'Order',  2, 'RemoveDC', 'on' );
                        case 2
                                %
                                % Perform polynomial detrending
                                %
                                fprintf('\n### Performing Polynomial Detrending...\n');
                                EEG  = nc_pop_polydetrend( EEG , 'Channels',  1:EEG.nbchan-1, 'Method', 'spline', ...
                                        'Windowsize',  3*EEG.srate, 'Windowstep',  1.5*EEG.srate ); % 3 seconds
                        otherwise
                                fprintf('\nWARNING: No high-pass filter applied.\n');
                end

                % #############################################################
                % RESAMPLE EEG DATA--------------------------------------------
                % #############################################################
                if params.downsample>1 && EEG.srate~=params.downsample
                        EEG = pop_resample( EEG, params.downsample );
                        EEG = eeg_checkset( EEG );
                end

                % #############################################################
                % LINE NOISE REMOVAL-------------------------------------------
                % #############################################################
                if removelinenoise % skip Cz
                        % Apply cleanLineNoise
                        fprintf('Removing line noise (50 or 60 Hz)...\n');
                        if kdata==1
                                lineNoiseIn = struct('lineNoiseMethod', 'clean', 'lineNoiseChannels', ...
                                        1:EEG.nbchan-1, 'Fs', EEG.srate, 'lineFrequencies', params.lineNoise.freq,...
                                        'p', 0.01, 'fScanBandWidth', 2, 'taperBandWidth', 2,...
                                        'taperWindowSize', 4, 'taperWindowStep', 4, ...
                                        'tau', 100, 'pad', 2, 'fPassBand', [0 EEG.srate/2], ...
                                        'maximumIterations', 10) ;
                        end
                        EEG = cleanLineNoise2(EEG, lineNoiseIn) ;
                end

                % #############################################################
                % LOW-PASS FILTER lpf_cutoff Hz (ERPLAB)-----------------------
                % #############################################################
                fprintf('Lowpass filtering the data...\n'); % skip Cz
                EEG  = pop_basicfilter( EEG,  1:EEG.nbchan-1 , 'Boundary', 'boundary',...
                        'Cutoff', lpf_cutoff,'Design', 'butter', 'Filter', 'lowpass',...
                        'Order',  2);

                % #############################################################
                % DETECT BAD CHANNELS HAPPE------------------------------------
                % #############################################################
                if params.badChans.rej
                        %origChans = EEG.chanlocs ;
                        [~, dataQC] = happe_detectBadChans(EEG, params, dataQC, ...
                                chanIDs, kdata) ;
                        lbl2int = regexp(dataQC{kdata,5},'E\d+\.?\d*','match');
                        badchanindx_happe = find(ismember({EEG.chanlocs.labels}, lbl2int));

                        % #######################################################
                        % DETECT BAD CHANNELS JAVIER-----------------------------
                        % #############################################################
                        fprintf('\n*** Searching again (Javier''s function) for bad channels...\n');
                        badchan     = nc_detectBadChannels(EEG, [1 0 0 0],[300 0 0 0], 25);
                        badchanindx_javier = find(badchan);

                        % merge bad channel indices
                        goodchanindx_all = 1:EEG.nbchan; % good channels (default all)
                        badchanindx_all  = unique([badchanindx_javier badchanindx_happe]) ;
                        lbch = length(badchanindx_all);

                        if lbch>0
                                % Javier's trick to carrier bad channel indices
                                EEG.comments = badchanindx_all;
                                fprintf('WARNING: Subject %s, has %g bad channels detected --> ', sname, lbch);
                                fprintf('Channel indices: [ ');
                                for kk=1:lbch
                                        fprintf('%g  ', badchanindx_all(kk));
                                end
                                fprintf(' ]\n');
                                % correct good channel indices
                                goodchanindx_all = goodchanindx_all(~ismember(goodchanindx_all, badchanindx_all ));

                                % correct bad chan data metric
                                bchanstr = sprintfc('E%d', badchanindx_all);
                                dataQC{kdata,3} = dataQC{kdata,2}-lbch;
                                dataQC{kdata,4} = dataQC{kdata,3}/dataQC{kdata,2}*100;
                                dataQC{kdata,5} = sprintf('%s ',bchanstr{:});
                        end
                else
                        dataQC{kdata,3} = size(chanIDs,2) ;
                        dataQC{kdata,4} = dataQC{kdata,3}/dataQC{kdata,2}*100;
                        dataQC{kdata,5} = 'NA' ;
                        badchanindx_happe = [];
                        goodchanindx_all  = 1:EEG.nbchan; % good channels (default all)
                end

                % #############################################################
                % WAVELET THRESHOLDING-----------------------------------------
                % #############################################################
                % wavelet threshold using happe_wavThresh's code
                EEG.data = double(EEG.data);
                %data2 = EEG.data;
                if EEG.srate > 500
                        wavLvl = 10;
                elseif EEG.srate > 250 && EEG.srate <= 500
                        wavLvl = 9;
                elseif EEG.srate <=250
                        wavLvl = 8;
                end
                wdata = reshape(EEG.data, size(EEG.data, 1), [])';
                ThresholdRule = 'Hard' ;
                artifacts = wdenoise(wdata, wavLvl, 'Wavelet', 'coif4', 'DenoisingMethod', ...
                        'Bayes', 'ThresholdRule', ThresholdRule, 'NoiseEstimate', ...
                        'LevelDependent')' ;

                preEEG      = reshape(EEG.data, size(EEG.data,1), []) ;
                postEEG     = preEEG - artifacts ;
                EEG.data    = postEEG ;
                EEG.setname = 'wavcleanedEEG' ;

                if ~isempty(goodchanindx_all)
                        % WAVELETING QC METRICS: Assesses the performance of wavelet thresholding.
                        % Only good channels
                        wavMeans = assessPipelineStep('wavelet thresholding', preEEG(goodchanindx_all,:), ...
                                postEEG(goodchanindx_all,:), wavMeans, EEG.srate, params.QCfreqs) ;
                        dataQC{kdata, 6} = var(postEEG(goodchanindx_all,:), 0, 'all')/var(preEEG(goodchanindx_all,:), ...
                                1, 'all')*100 ;
                else
                        dataQC{kdata, 6} = 'too bad';
                end

                % #############################################################
                % INTERPOLATE BAD CHANNELS----(using cleaner data)-----------
                % #############################################################
                if params.segment.interp
                        if length(badchanindx_all)<(EEG.nbchan*propbadch) && ~isempty(badchanindx_all)
                                EEG = eeg_interp(EEG, badchanindx_all, 'spherical');
                                % Chans Interpolated
                                dataQC{kdata,7} = length(badchanindx_all) ;
                        else
                                % Chans Interpolated
                                dataQC{kdata,7} = ['Too_many (' num2str(length(badchanindx_all)) ')' ] ;
                        end
                else
                        %Chans Interpolated
                        dataQC{kdata,7} = 'NA' ;
                end
                % #############################################################
                % AVERAGE RE-REFERENCE
                % #############################################################
                if params.reref.on
                        fprintf('Re-referencing...\n') ;
                        EEG = pop_reref(EEG, [], 'keepref', 'on') ;
                        refstr = 'avgreref';
                else
                        refstr = 'noreref';
                end

                % eegplot( EEG.data, 'srate', EEG.srate, 'title', 'Black = channel before denoising; red = after denoising -- eegplot()', ...
                %           'limits', [EEG.xmin EEG.xmax]*1000, 'data2', data2);

                %###############################################################
                % COUNT OCCURENCE OF EVENT CODES CT02 AND IT01 AND FIND SPLITTING POINT
                %###############################################################
                totalpnts = EEG.pnts;
                evType = categorical({EEG.event.type});
                [Ntype, LabelType] = histcounts(evType);
                ct02_count = Ntype(string(LabelType)=="ct02") ;
                it01_count = Ntype(string(LabelType)=="it01") ;

                % remove all event codes but ct02 and it01
                evmask    = string({EEG.event.type})=="it01" | string({EEG.event.type})=="ct02";
                Auxevent  = EEG.event;
                isolatedEvent = EEG.event(evmask);

                % block counter (count 12 --> 12*4 = 48 even codes)
                evType = categorical({isolatedEvent.type});
                iev = 2;
                counblock = 1;
                while iev<=length(evType) && counblock<=12

                        if evType(iev)~=evType(iev-1)
                                counblock = counblock+1;
                        end
                        iev = iev+1;
                end
                % split point calculation (to use right before saving the datasets below)
                point2split = isolatedEvent(iev-1).latency - 200*EEG.srate/1000; % preserve -200 ms
                %###############################################################

                % #############################################################
                % SAVE PROCESSED CONTINUOUS DATASET --- PRE portion ------------------------
                % #############################################################

                EEGaux = EEG; % keep the whole dataste for later

                EEG = pop_select( EEG, 'point',[0 point2split] ); % pre portion
                EEG = eeg_checkset( EEG );
                fprintf('\n%s\n\n', repmat('#',1, 120));
                fprintf('\nSaving the "PRE" portion of the current dataset...Num Chan = %g\n', EEG.nbchan);
                fprintf('\n%s\n\n', repmat('#',1, 120));
                if true_CB == 1
                        CB_str = 'CB1';
                elseif true_CB == 2
                        CB_str = 'CB2';
                else
                        CB_str = 'CB_UNKNOWN';
                end

                %
                % NOTE: So, now the filename should contain all the relevant
                % info for bein used by scripts 2 and 3. In this way, the spreadsheet with
                % the subjects' info will not be requiered by scripts 2 and 3.
                % JLC
                %
                session = 'PRE';
                % folder name to write the final processed dataset
                pathname_write = fullfile(path_source, 'PREPROC_ERPTRAIN_01', upper(session), [TaskName '_' subj_ID]);

                if exist(pathname_write, 'dir')==0
                        % if folder doesn't exist then make it
                        mkdir(pathname_write);
                end
                EEG.setname = sprintf('%s_%s_%s_%s_%s_%s_Preproc_ERPTRAIN_Wavclean_%s',upper(TaskName), upper(session), subj_ID, subj_AGE, subj_observa, CB_str, refstr);
                EEG = pop_saveset( EEG,  'filename', [EEG.setname  '.set'], 'filepath', pathname_write);

                clear EEG % clear current splitted portion

                % #############################################################
                % SAVE PROCESSED CONTINUOUS DATASET --- POST portion ------------------------
                % #############################################################

                EEG = EEGaux; % recover the whole processed EEG

                EEG = pop_select( EEG, 'point',[point2split EEG.pnts] ); % post portion
                EEG = eeg_checkset( EEG );
                fprintf('\n%s\n\n', repmat('#',1, 120));
                fprintf('\nSaving the "POST" portion of the current dataset...Num Chan = %g\n', EEG.nbchan);
                fprintf('\n%s\n\n', repmat('#',1, 120));

                %
                % NOTE: So, now the filename should contain all the relevant
                % info for bein used by scripts 2 and 3. In this way, the spreadsheet with
                % the subjects' info will not be requiered by scripts 2 and 3.
                % JLC
                %
                session = 'POST';
                % folder name to write the final processed dataset
                pathname_write = fullfile(path_source, 'PREPROC_ERPTRAIN_01', upper(session), [TaskName '_' subj_ID]);

                if exist(pathname_write, 'dir')==0
                        % if folder doesn't exist then make it
                        mkdir(pathname_write);
                end
                EEG.setname = sprintf('%s_%s_%s_%s_%s_%s_Preproc_ERPTRAIN_Wavclean_%s',upper(TaskName), upper(session), subj_ID, subj_AGE, subj_observa, CB_str, refstr);
                EEG = pop_saveset( EEG,  'filename', [EEG.setname  '.set'], 'filepath', pathname_write);

                clear EEG EEGaux

                eval_str = 'PASSED';
%         catch ME2
%                 eval_str =  ME2.message;
%                 fprintf('\n%s\n\n', repmat('9',1, 120));
%                 fprintf('\n\nOops!: something went wrong....[ %s ]\n\n', eval_str) ;
%                 fprintf('\n%s\n\n', repmat('9',1, 120));
%                 beep
%                 [dataQC{kdata,1:7}] = deal('crashed  ');
% 
%         end
        dataQC{kdata,8}  = sname; % filename
        dataQC{kdata,9}  = Age_from_file;
        dataQC{kdata,10} = true_AGE;
        dataQC{kdata,11} = Obs_from_file;
        dataQC{kdata,12} = true_Obs;
        dataQC{kdata,13} = CB_from_file;
        dataQC{kdata,14} = true_CB ;
        dataQC{kdata,15} = ct02_count;
        dataQC{kdata,16} = it01_count;
        dataQC{kdata,17} = totalpnts; % EEG.pnts
        dataQC{kdata,18} = point2split; % sample position for splitting
        dataQC{kdata,19} = eval_str ;
end
disp('Processing Done!')
fprintf('\nNow saving the table...\n') ;
T = cell2table(dataQC,'VariableNames',dataQCnames)
% reorganize variables of the table
T = movevars(T,{'Filename','Age_from_file','True_age','Obs_from_file',...
        'True_Obs','CB_from_file','True_CB','ct02_count',...
        'it01_count','Amount_Samples','Splitting_sample','QC'},'Before','File_Length_in_Seconds');
% save the table
taildate = sprintf('%12.0f',now*10^6);
tablename = fullfile(path_source, sprintf('OUTPUT_ERPTRAIN_SCRIPT_1_%s_%s_Preproc_Wavclean_%s_%s.xlsx',upper(TaskName), upper(session), refstr, taildate));
writetable(T, tablename);
save(fullfile(paramFilePath, sprintf('%s_%s.mat',paramFileName, taildate)), 'params');
disp('All Done!')
dt2 = datetime;
fprintf('\n* Script --> start: %s , end: %s , duration: %s\n',dt1,dt2,dt2-dt1);
