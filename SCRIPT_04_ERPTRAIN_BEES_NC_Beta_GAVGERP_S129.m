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
% close all
dt1 = datetime;
% ##########################################################################
% ##########################  EDIT STARTS  #################################
% ##########################################################################


TaskName = 'BEES';
session   = {'pre','post'}; % 'post' or/and 'pre'

% folder where you have the raw data
path_source   = '/Volumes/WORKBENCH/LISADATA/BEES_ERR_mff files';

workingchannels = 129; % correct number of (final) channels
workingsrate    = 500; % correct sample rate

% ##########################################################################
% ##########################   EDIT ENDS  ##################################
% ##########################################################################

jsession6  =  1;
jsession9  =  1;
jsession12 =  1;

for isession = 1:length(session) % in case more than

        currsession = session{isession};

        % PREPARE DATASET LIST TO BE PROCESSED-------------------------------------
        pathname_read  = fullfile(path_source, 'ERP_ERPTRAIN_02', upper(currsession)); % ### EDIT THIS ###
        pathname_write = fullfile(path_source, 'GAVGERP_ERPTRAIN', upper(currsession));% ### EDIT THIS ###
        if exist(pathname_write, 'dir')==0
                % if folder doesn't exist then make it
                mkdir(pathname_write);
        end
        a = dir(pathname_read);
        b = regexpi({a.name}','BEES_\d+','match'); % ### EDIT THIS ###
        folder  = [b{:}];
        nfolder = length(folder);

        hdata6 = 1;
        hdata9 = 1;
        hdata12 = 1;

        for kfolder=1:nfolder

                pathname_read_local  = fullfile(pathname_read, folder{kfolder});

                fprintf('\n---------------------------------------------\n');
                fprintf('Searching for existing datasets...\n');
                % c = regexpi(b,'BEES\s*[\w+|\d+].*(?=\.)','match');
                a = dir(pathname_read_local);
                b = regexpi({a.name}','BEES\s*[\w+|\d+].*(?=\.erp)','match'); % *.erp % ### EDIT THIS ###
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
                                                subj_AGE  = str2num(subj_info{2}); % age from filename
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
                                ERP = pop_loaderp( 'filename', [sname '.erp'], 'filepath', pathname_read_local);

                                % CHECK CHANNEL NUMBER (CHCK 1)--------------------------------------
                                if ERP.nchan~=workingchannels
                                        error('wrong channel number')
                                end

                                % Art det
                                MPD = getardetection(ERP, 1);

                                if MPD<=50 && isnumeric(subj_AGE) % only ERPs having better than 50% artifact detection will be grand averaged

                                        % create new bin to get CT-IT
                                        ERP = pop_binoperator( ERP, {  'b3 = b2-b1 label diffwave'}); % ### EDIT THIS ###

                                        %
                                        % Add new channels (ROIs)
                                        %
                                        % % Regions of interest
                                        % LOT = [58,59,64];
                                        % LO  = [65,69,70];
                                        % MO  = [74,75,82];
                                        % RO  = [83,89,90];
                                        % ROT = [91,95,96];

                                        chanFormulas = {...
                                                'ch107 = (ch58 + ch59 + ch64)/3 label = LOT',...
                                                'ch108 = (ch65 + ch69 + ch70)/3 label = LO',...
                                                'ch109 = (ch74 + ch75 + ch82)/3 label = MO',...
                                                'ch110 = (ch83 + ch89 + ch90)/3 label = RO',...
                                                'ch111 = (ch91 + ch95 + ch96)/3 label = ROT'...
                                                }; % ### EDIT THIS ###

                                        ERP = pop_erpchanoperator( ERP, chanFormulas, 'KeepLocations',  1, 'Warning', 'off' );

                                        switch subj_AGE
                                                case 6
                                                        ALLERP6(hdata6) = ERP;
                                                        hdata6 = hdata6 + 1; % pointer of good ERPs
                                                case 9
                                                        ALLERP9(hdata9) = ERP;
                                                        hdata9 = hdata9 + 1; % pointer of good ERPs
                                                case 12
                                                        ALLERP12(hdata12) = ERP;
                                                        hdata12 = hdata12 + 1; % pointer of good ERPs
                                                otherwise
                                                        fprintf('\n###%s was excluded from the Grand Average...\n', [sname '.erp']);

                                        end
                                else
                                        fprintf('\n###%s was excluded from the Grand Average...\n', [sname '.erp']);
                                end

                        catch ME2
                                eval_str =  ME2.message;
                                fprintf('\n%s\n\n', repmat('9',1, 120));
                                fprintf('\n\nOops!: something went wrong....[ %s ]\n\n', eval_str) ;
                                fprintf('\n%s\n\n', repmat('9',1, 120));
                                beep
                        end

                        clear ERP
                end
        end

        if hdata6>1 % only ERPs 6-month having better than 50% artifact detection will be grand averaged

                %
                % Get Grand Averages (weighted average)
                %
                GAVGERP6(jsession6) = pop_gaverager( ALLERP6 , 'Criterion',50, 'Erpsets', 1:length(ALLERP6), 'ExcludeNullBin', 'on', 'SEM', 'on' );

                %
                % Save Grand Average(s)
                %
                fprintf('\n###Saving Grand Average 6 month..\n');

                erpname6 = sprintf('GAVG_6m_%s_%s_ERPTRAIN',upper(TaskName), upper(currsession));
                GAVGERP6(jsession6).erpname = erpname6;

                pop_savemyerp(GAVGERP6(jsession6), 'erpname', erpname6, 'filename', [erpname6 '.erp'], 'filepath', path_source);
                fprintf('Processed erpset was saved!\n\n');

                clear folder ALLERP6
                jsession6 = jsession6 + 1; % pointer of good session
        end
        if hdata9>1 % only ERPs 9-month having better than 50% artifact detection will be grand averaged

                %
                % Get Grand Averages (weighted average)
                %
                GAVGERP9(jsession9) = pop_gaverager( ALLERP9 , 'Criterion',50, 'Erpsets', 1:length(ALLERP9), 'ExcludeNullBin', 'on', 'SEM', 'on' );

                %
                % Save Grand Average(s)
                %
                fprintf('\n###Saving Grand Average 9 month..\n');

                erpname9 = sprintf('GAVG_9m_%s_%s_ERPTRAIN',upper(TaskName), upper(currsession));
                GAVGERP9(jsession9).erpname = erpname9;

                pop_savemyerp(GAVGERP9(jsession9), 'erpname', erpname9, 'filename', [erpname9 '.erp'], 'filepath', path_source);
                fprintf('Processed erpset was saved!\n\n');

                clear folder ALLERP9
                jsession9 = jsession9 + 1; % pointer of good session
        end
        if hdata12>1 % only ERPs 12-month having better than 50% artifact detection will be grand averaged

                %
                % Get Grand Averages (weighted average)
                %
                GAVGERP12(jsession12) = pop_gaverager( ALLERP12 , 'Criterion',50, 'Erpsets', 1:length(ALLERP12), 'ExcludeNullBin', 'on', 'SEM', 'on' );

                %
                % Save Grand Average(s)
                %
                fprintf('\n###Saving Grand Average 12 month..\n');

                erpname12 = sprintf('GAVG_12m_%s_%s_ERPTRAIN',upper(TaskName), upper(currsession));
                GAVGERP12(jsession12).erpname = erpname12;

                pop_savemyerp(GAVGERP12(jsession12), 'erpname', erpname12, 'filename', [erpname12 '.erp'], 'filepath', path_source);
                fprintf('Processed erpset was saved!\n\n');

                clear folder ALLERP12
                jsession12 = jsession12 + 1; % pointer of good session
        end
end

if jsession6>1
        %
        % Append GAVG 6 months ERPs coming from all specified sessions
        %
        APPGAVGERP6 = pop_appenderp( GAVGERP6 , 'Erpsets', 1:length(GAVGERP6), 'Prefixes', 'erpname' );

        %
        % Get appended 6 months Grand Averages if more tham one session was specified
        %
        erpname6 = sprintf('APPE_GAVG_6m_%s_ERPTRAIN_ALL',upper(TaskName));
        pop_savemyerp(APPGAVGERP6, 'erpname', erpname6, 'filename', [erpname6 '.erp'], 'filepath', path_source);
        fprintf('Appended GAVG 6 months was saved!\n\n');
end

if jsession9>1
        %
        % Append GAVG 9 months ERPs coming from all specified sessions
        %
        APPGAVGERP9 = pop_appenderp( GAVGERP9 , 'Erpsets', 1:length(GAVGERP9), 'Prefixes', 'erpname' );

        %
        % Get appended 9 months Grand Averages if more tham one session was specified
        %
        erpname9 = sprintf('APPE_GAVG_9m_%s_ERPTRAIN_ALL',upper(TaskName));
        pop_savemyerp(APPGAVGERP9, 'erpname', erpname9, 'filename', [erpname9 '.erp'], 'filepath', path_source);
        fprintf('Appended GAVG 9 months was saved!\n\n');
end

if jsession12>1
        %
        % Append GAVG 12 months ERPs coming from all specified sessions
        %
        APPGAVGERP12 = pop_appenderp( GAVGERP12 , 'Erpsets', 1:length(GAVGERP12), 'Prefixes', 'erpname' );

        %
        % Get appended 12 months Grand Averages if more tham one session was specified
        %
        erpname12 = sprintf('APPE_GAVG_12m_%s_ERPTRAIN_ALL',upper(TaskName));
        pop_savemyerp(APPGAVGERP12, 'erpname', erpname12, 'filename', [erpname12 '.erp'], 'filepath', path_source);
        fprintf('Appended GAVG 12 months was saved!\n\n');
end


if jsession6>1 && jsession9>1 && jsession12>1
        %
        % Append ALL GAVG ERPs (6,9, and 12 months together)
        %
        APPGAVGERPALL(1)= APPGAVGERP6;
        APPGAVGERPALL(2)= APPGAVGERP9;
        APPGAVGERPALL(3)= APPGAVGERP12;

        ERP = pop_appenderp( APPGAVGERPALL , 'Erpsets', 1:3);

        %
        % Get appended all gavgs
        %
        erpname = sprintf('APPE_GAVG_ALLAGES_%s_ERPTRAIN_ALL',upper(TaskName));
        pop_savemyerp(ERP, 'erpname', erpname, 'filename', [erpname '.erp'], 'filepath', path_source);
        fprintf('Appended GAVG all ages was saved!\n\n');

else
        fprintf('ERROR: Final appending of GAVGs did not work ...\n\n');
end










