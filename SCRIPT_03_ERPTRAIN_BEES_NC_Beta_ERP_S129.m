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
session  = 'post'; % 'post' or 'pre'

% folder where you have the raw data
path_source   = '/Volumes/WORKBENCH/LISADATA/BEES_ERR_mff files';

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

% PREPARE DATASET LIST TO BE PROCESSED-------------------------------------
pathname_read  = fullfile(path_source, 'EPOCHED_ERPTRAIN_02', upper(session));
pathname_write = fullfile(path_source, 'ERP_ERPTRAIN_02', upper(session));
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

                        % get good epoch indices per bin to attach this info to the
                        % ERP structure
                        goodepochsAll = bitand(not(sum(EEG.reject.rejmanualE,1)>0), not(EEG.reject.rejmanual));
                        nbin = EEG.EVENTLIST.nbin;
                        FinalSelectedEpochs = zeros(nbin, EEG.trials);
                        for ibin = 1:nbin
                                FinalSelectedEpochs(ibin, epoch4bin(EEG, ibin)) = 1;
                                FinalSelectedEpochs(ibin, :) = bitand(FinalSelectedEpochs(ibin, :), goodepochsAll);
                        end

                        % Art det
                        [~, MPD] = getardetection(EEG, 1);
                        S(hdata).name = sname;
                        S(hdata).reject = MPD;

                        %
                        % Get averaged ERP
                        %
                        ERP = pop_averager( EEG , 'Criterion', 'good', 'ExcludeBoundary', 'on', 'SEM', 'on' );
                        % attach good epoch indices per bin
                        ERP.ntrials.FinalSelectedEpochs = FinalSelectedEpochs;

                        S(hdata).ntrialaccepted = ERP.ntrials.accepted;
                        S(hdata).ntrialrejected = ERP.ntrials.rejected;

                        %
                        % Save averaged ERP
                        %
                        fprintf('\n###Saving erpset...\n');
                        fprintf('\nSaving subject %s \n', sname);

                        if subj_CB == 1
                                CB_str = 'CB1';
                        elseif subj_CB == 2
                                CB_str = 'CB2';
                        else
                                CB_str = 'CB_UNKNOWN';
                        end

                        erpname = sprintf('%s_%s_%s_%s_%s_%s_Preproc_ERPTRAIN_ERP',upper(TaskName), upper(session), subj_ID, subj_AGE, subj_observa, CB_str);
                        ERP     = pop_savemyerp(ERP, 'erpname', erpname, 'filename', [erpname '.erp'], 'filepath', pathname_write_local);
                        fprintf('Processed erpset was saved!\n\n');
                        eval_str = 'PASSED';
                catch ME2
                        eval_str =  ME2.message;
                        fprintf('\n%s\n\n', repmat('9',1, 120));
                        fprintf('\n\nOops!: something went wrong....[ %s ]\n\n', eval_str) ;
                        fprintf('\n%s\n\n', repmat('9',1, 120));
                        beep

                        S(hdata).name = sname;
                        S(hdata).reject = "crashed";
                        S(hdata).ntrialaccepted = "crashed";
                        S(hdata).ntrialrejected = "crashed";
                end

                hdata = hdata + 1;
                clear ERP EEG
        end
end
disp('Done!')
T = struct2table(S);
% reorganize variables of the table
taildate = sprintf('%12.0f',now*10^6);
tablename = fullfile(path_source, sprintf('GrandTable_%s_%s_ERPTRAIN_ERP_%s.xlsx',upper(TaskName), upper(session), taildate));
writetable(T, tablename);
dt2 = datetime;
fprintf('\n* Script --> start: %s , end: %s , duration: %s\n',dt1,dt2,dt2-dt1);