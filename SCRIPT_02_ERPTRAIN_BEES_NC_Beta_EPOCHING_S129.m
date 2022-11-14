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
TaskName = 'BEES';
session   = 'pre'; % 'post' or 'pre'

% folder where you have the raw data
path_source   = '/Volumes/WORKBENCH/LISADATA/BEES_ERR_mff files';

% Optional:
% folder where you want to save individual conditions as a mat file
saveIndividualMatFiles = true; % true or false
path_matfiles = '/Users/Javier/Documents/MATLAB/LISA_DATA/MATFILES';

% For event code re-mapping (creates bins without binlister)
eventmappath = '/Volumes/WORKBENCH/LISADATA/BCDLAB_FTAG_beta001/Infants/Newencode';
eventmapfile = 'event_remapping_erptrain.txt';

% epoching window in msec
epochwin = [-200  800];

workingchannels = 129; % correct number of (final) channels
workingsrate    = 500; % correct sample rate

%   8 = 'E8'
%   9 = 'E9'
%  14 = 'E14'
%  15 = 'E15'
%  20 = 'E21'
%  21 = 'E22'
%  24 = 'E25'
eblnkchan = [8 9 14 15 18 20 21 24];  % eye movement channels

%
% if you want this script to set your Matlab path -automatically- by adding the required
% packages then make loadPackagesIntoMatlab = true
%
loadPackagesIntoMatlab = false;

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
% -------------------------------------------------------------------------
if saveIndividualMatFiles
        % check existency of folder to save individual mat files (if needed)
        if exist(path_matfiles, 'dir')==0
                % if folder doesn't exist then make it
                mkdir(path_matfiles);
        end
end

% PREPARE DATASET LIST TO BE PROCESSED-------------------------------------
pathname_read  = fullfile(path_source, 'PREPROC_ERPTRAIN_01', upper(session));
pathname_write = fullfile(path_source, 'EPOCHED_ERPTRAIN_02', upper(session));
if exist(pathname_write, 'dir')==0
        % if folder doesn't exist then make it
        mkdir(pathname_write);
end
a = dir(pathname_read);
b = regexpi({a.name}','BEES_\d+','match');
folder  = [b{:}];
nfolder = length(folder);
hdata = 1;

for kfolder=1:nfolder

        pathname_read_local  = fullfile(pathname_read, folder{kfolder});
        % folder name to write the epoched dataset
        pathname_write_local = fullfile(pathname_write, folder{kfolder});
        if exist(pathname_write_local, 'dir')==0
                % if folder doesn't exist then make it
                mkdir(pathname_write_local);
        end
        fprintf('\n---------------------------------------------\n');
        fprintf('Searching for existing datasets...\n');
        % c = regexpi(b,'BEES\s*[\w+|\d+].*(?=\.)','match');
        a = dir(pathname_read_local);
        b = regexpi({a.name}','BEES\s*[\w+|\d+].*(?=\.set)','match');
        if isempty(b)
                nfile = 0;
        else
                fileArray  = [b{:}];
                nfile = length(fileArray);
                if kfolder==1
                        % INITIALIZE QUALITY REPORT METRICS----------------------------------------
                        fprintf('Initializing report metrics...\n') ;
                        % DATA QUALITY METRICS: create a variable holding all the names of each
                        % metric and a variable to hold the metrics for each file.
                        dataQCnames = {'Epoch_Length_in_mSec', 'Number_User-Selected_Chans', ...
                                'Number_Segs_Pre-Seg_Rej', 'Number_Segs_Post-Seg_Rej', ...
                                'Percent_Segs_Post-Seg_Rej'} ;
                        dataQC     = cell(nfile, length(dataQCnames)) ;
                end
        end

        fprintf('\nProcessing dataset %g of %g ... Be patient please...\n\n', kfolder, nfolder);

        % DATASET LOOP-------------------------------------------------------------
        for kdata=1:nfile
                try
                        sname   = fileArray{kdata}; % current filename

                        % CB from filename
                        get_CB = contains(sname,'cb2','IgnoreCase',true);
                        if get_CB
                                subj_CB = 2; % trust it this time
                        else
                                subj_CB = 1; % trust it this time
                        end

                        % subject ID from filename.  Trust it this time
                        subj_info = regexpi(sname,'_(\d{3}+)_(\d{1}+)*[_]*(\d{1}+)*','tokens'); % get ID and age and observation (if any)
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

                        % LOAD PREPROC DATA------------------------------------------------------
                        fprintf('\n%s\n\n', repmat('#',1, 150));
                        fprintf('Loading %s.set... Please wait\n', sname);
                        fprintf('\n%s\n\n', repmat('#',1, 150));
                        EEG = pop_loadset('filename',[sname '.set'],'filepath', pathname_read_local);

                        % CHECK CHANNEL NUMBER (CHCK 1)--------------------------------------
                        if EEG.nbchan~=workingchannels
                                error('wrong channel number')
                        end

                        % channels over scalp (not eblnkchan)
                        scalpchan = find(~ismember(1:EEG.nbchan, eblnkchan));

                        % CREATE EVENTLIST AND BINS DIRECTLY, BY MEANS OF CODE REMAPPING (IN ONE STEP, NO BINLISTER)
                        % B1(it01) B2(iu02) B3(ct03) B4(cu04) B5(un05)
                        EEG  = pop_editeventlist( EEG , 'BoundaryNumeric', { -99}, 'BoundaryString', { 'boundary' },...
                                'List', fullfile(eventmappath,eventmapfile), 'SendEL2', 'EEG', 'UpdateEEG', 'code', 'Warning', 'off' );

                        EEG = pop_epochbin( EEG , epochwin,  [0 100]);

                        % Report file length in msec'
                        dataQC{hdata, 1} = 1000*(EEG.xmax-EEG.xmin) ;
                        dataQC{hdata, 2} = workingchannels;

                        %
                        % Artifact detection: 1st round
                        %
                        EEG  = pop_artextval( EEG , 'Channel',  scalpchan, 'Flag',  1, 'Threshold', [ -120 120],...
                                'Twindow', [ -200 200] );

                        % ###########################################################
                        % check for remaining bad channels that takes most of
                        % the artifact detection (usually 1 damn bad channel)
                        % However, this time pick only 3 channels that explain
                        % more than 50% of the artifacts

                        Frej = fieldnames(EEG.reject);
                        sfields1 = regexpi(Frej, '\w*E$', 'match');
                        sfields2 = [sfields1{:}];
                        nfield   = length(sfields2);
                        histE    = zeros(EEG.nbchan, EEG.trials);
                        for i=1:nfield
                                fieldnameE = char(sfields2{i});
                                if ~isempty(EEG.reject.(fieldnameE))
                                        histE = histE | [EEG.reject.(fieldnameE)]; %electrodes
                                end
                        end
                        histeEF = sum(histE,2)';
                        Ttop = EEG.trials;
                        badchanindx = find([(histeEF/Ttop)*100]>50);
                        if ~isempty(badchanindx)
                                if length(badchanindx)>5
                                        badchanindx = badchanindx(1:5);
                                end
                                fprintf('\n### Interpolating %g channels...\n', length(badchanindx));
                                EEG = eeg_interp(EEG, badchanindx, 'spherical');
                        end
                        % ###########################################################

                        %
                        % Artifact detection: 2nd round
                        %
                        EEG  = pop_artextval( EEG , 'Channel',  scalpchan, 'Flag',  1, 'Threshold', [ -120 120],...
                                'Twindow', [ -200 200] );

                        % Art det
                        [~, tprej, acce, rej] = pop_summary_AR_eeg_detection(EEG, 'none');
                        dataQC{hdata, 3} = EEG.trials;
                        dataQC{hdata, 4} = sum(acce);
                        dataQC{hdata, 5} = tprej;

                        if hdata==1
                                bindesc = {EEG.EVENTLIST.bdf.description};
                        end
                        accrej = [];
                        for ibin=1:length(bindesc)
                                if kfolder==1 && kdata==1
                                        dataQCnames = [ dataQCnames ['ACCE_' bindesc{ibin}]];
                                        dataQCnames = [ dataQCnames ['REJE_' bindesc{ibin}]];
                                end
                                dataQC{hdata, 6+2*(ibin-1)} = acce(ibin) ;
                                dataQC{hdata, 7+2*(ibin-1)} = rej(ibin);
                        end

                        % ###########################################################
                        %
                        % Save binepoched  dataset
                        %
                        % ###########################################################

                        if subj_CB == 1
                                CB_str = 'CB1';
                        elseif subj_CB == 2
                                CB_str = 'CB2';
                        else
                                CB_str = 'CB_UNKNOWN';
                        end
                        setname_set = sprintf('%s_%s_%s_%s_%s_%s_ERPTRAIN_EPOCHED',upper(TaskName), upper(session), subj_ID, subj_AGE, subj_observa, CB_str);

                        if saveIndividualMatFiles
                                % ###########################################################
                                %
                                % Export separated .mat files in sake of data democracy
                                %
                                % Epochs corresponding to each individual event code are
                                % saved as separated -mat files
                                % ###########################################################
                                ibin = 1:EEG.EVENTLIST.nbin;
                                for pbin=1:length(ibin)
                                        % get epoch indices from "good" trials for a specific bin
                                        EpochIndex  = getepochindex6(EEG, 'Bin', ibin(pbin), 'Nepoch', 'amap','Artifact','good');
                                        EpochIndex  = cell2mat(EpochIndex);
                                        if isempty(EpochIndex)
                                                fprintf('Ooops! No good epochs were found for %s...\n', EEG.EVENTLIST.bdf(pbin).description);
                                        else
                                                fprintf('Saving separated mat files for %s having %g epochs...\n', EEG.EVENTLIST.bdf(pbin).description, length(EpochIndex));
                                                epochData   = EEG.data(:,:, EpochIndex);
                                                setname_mat = sprintf('%s_%s', setname_set, EEG.EVENTLIST.bdf(pbin).description);
                                                save(fullfile( path_matfiles,[setname_mat '.mat']),'epochData', '-v7.3');
                                                clear epochData
                                        end
                                end
                        end
                        fprintf('Saving dataset...Num Chan = %g\n', EEG.nbchan);
                        EEG.setname = setname_set;
                        EEG = pop_saveset( EEG,  'filename', [EEG.setname  '.set'], 'filepath', pathname_write_local);

                        clear EEG
                        eval_str = 'PASSED';
                catch ME2
                        eval_str =  ME2.message;
                        fprintf('\n%s\n\n', repmat('9',1, 120));
                        fprintf('\n\nOops!: something went wrong....[ %s ]\n\n', eval_str) ;
                        fprintf('\n%s\n\n', repmat('9',1, 120));
                        beep
                        [dataQC{hdata,1:7+2*(length(bindesc)-1)}] = deal('crashed  ');
                end
                dataQC{hdata,7+2*(length(bindesc)-1)+1} = sname ;        % filename
                dataQC{hdata,7+2*(length(bindesc)-1)+2} = num2str(badchanindx) ; % bad ch
                dataQC{hdata,7+2*(length(bindesc)-1)+3} = subj_ID ;      % ID
                dataQC{hdata,7+2*(length(bindesc)-1)+4} = subj_AGE ;     % Age
                dataQC{hdata,7+2*(length(bindesc)-1)+5} = subj_observa ; % Observa
                dataQC{hdata,7+2*(length(bindesc)-1)+6} = eval_str ;     % QC

                hdata = hdata + 1;
        end
end
dataQCnames = [ dataQCnames 'Filename' 'Chans_Interpol' 'ID' 'AGE' 'OBSERVA' 'QC'];
disp('Done!')
T = cell2table(dataQC,'VariableNames',dataQCnames);
% reorganize variables of the table
T = movevars(T,{'Filename' 'ID' 'AGE' 'OBSERVA' 'QC'},'Before','Epoch_Length_in_mSec');
T = movevars(T,{'Chans_Interpol'},'After','Number_User-Selected_Chans')
taildate = sprintf('%12.0f',now*10^6);
tablename = fullfile(path_source, sprintf('GrandTable_%s_%s_ERPTRAIN_EPOCHED_%s.xlsx',upper(TaskName), upper(session), taildate));
writetable(T, tablename);
dt2 = datetime;
dt2 = datetime;
fprintf('\n* Script --> start: %s , end: %s , duration: %s\n',dt1,dt2,dt2-dt1);