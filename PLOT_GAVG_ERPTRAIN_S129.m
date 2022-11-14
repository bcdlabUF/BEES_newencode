% Load GAVG ERPs, plot IT, CT, and diff waves bins
% (ERPLAB installed is required)
%
% ERP-TRAIN task
% This are the bins withing the grand-average
%
% BIN1 = GAVG_6m_BEES_PRE_ERPTRAIN : it01
% BIN2 = GAVG_6m_BEES_PRE_ERPTRAIN : ct02
% BIN3 = GAVG_6m_BEES_PRE_ERPTRAIN : diffwave
% BIN4 = GAVG_6m_BEES_POST_ERPTRAIN : it01
% BIN5 = GAVG_6m_BEES_POST_ERPTRAIN : ct02
% BIN6 = GAVG_6m_BEES_POST_ERPTRAIN : diffwave
% BIN7 = GAVG_9m_BEES_PRE_ERPTRAIN : it01
% BIN8 = GAVG_9m_BEES_PRE_ERPTRAIN : ct02
% BIN9 = GAVG_9m_BEES_PRE_ERPTRAIN : diffwave
% BIN10 = GAVG_9m_BEES_POST_ERPTRAIN : it01
% BIN11 = GAVG_9m_BEES_POST_ERPTRAIN : ct02
% BIN12 = GAVG_9m_BEES_POST_ERPTRAIN : diffwave
% BIN13 = GAVG_12m_BEES_PRE_ERPTRAIN : it01
% BIN14 = GAVG_12m_BEES_PRE_ERPTRAIN : ct02
% BIN15 = GAVG_12m_BEES_PRE_ERPTRAIN : diffwave
% BIN16 = GAVG_12m_BEES_POST_ERPTRAIN : it01
% BIN17 = GAVG_12m_BEES_POST_ERPTRAIN : ct02
% BIN18 = GAVG_12m_BEES_POST_ERPTRAIN : diffwave
%
%
% This are the channels containing the ROIs
% 107 = LOT
% 108 = LO
% 109 = MO
% 110 = RO
% 111 = ROT
%
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

clc
loadgavgerpOp = false;
if exist('ERP','var')
    answer = questdlg('What would you like to plot?', ...
        'Hello','a new GAVG ERP file', 'just continue working','Cancel','Cancel');
    switch answer
        case 'a new GAVG ERP file'
            disp('Coming right up!')
            loadgavgerpOp = true;
        case 'just continue working'
            disp('You betcha!')
        case 'Cancel'
            disp('User selected Cancel');
            return
    end
else
    loadgavgerpOp = true;
end
if loadgavgerpOp
    % load ERP (all conditions) 
    [filegavgerp, pathgavgerp] = uigetfile('*.erp','GAVG ERP File Selector');
    if isequal(filegavgerp,0)
        disp('User selected Cancel');
        return
    end
    ERP = pop_loaderp( 'filename', filegavgerp, 'filepath', pathgavgerp);
end

ChanArray      = 107:111; 
BinArray_IT    = 1:3:16;
BinArray_CT    = 2:3:17;
BinArray_DIFF  = 3:3:18;
LineWidth = 2;
xscale    = [-200.0 798.0  -200:200:600]; % [xmin xmax xticks]
yscale    = [ -8.0  15.0   -8:2:16 ] ;    % [ymin ymax yticks]
YDir      = 'normal'; % 'normal' , 'reverse'
FontSizeChan  = 18; 
FontSizeLeg   = 14;
FontSizeTicks = 14;
baseline      = '-200 0';

% 
% ERP = pop_ploterps( ERP,  1:4,  107:111 , 'AutoYlim', 'on', 'Axsize', [ 0.05 0.08], 'BinNum', 'on', 'Blc', 'pre', 'Box', [ 5 1], 'ChLabel',...
%  'on', 'FontSizeChan',  12, 'FontSizeLeg',  12, 'FontSizeTicks',  12, 'LegPos', 'right', 'Linespec', {'k-' , 'r-' , 'b-' , 'g-' }, 'LineWidth',...
%   1, 'Maximize', 'on', 'Position', [ 103.714 29.6429 106.857 31.9286], 'Style', 'Classic', 'Tag', 'ERP_figure', 'Transparency',  0, 'xscale',...
%  [ -200.0 798.0   -200:200:600 ], 'YDir', 'normal' );


%
% plot ITs
%
ERP = pop_ploterps( ERP,  BinArray_IT ,  ChanArray , 'Axsize', [ 0.05 0.08], 'BinNum', 'on', 'Blc', baseline, 'Box', [ 5 1], 'ChLabel', 'on',...
 'FontSizeChan',  FontSizeChan, 'FontSizeLeg',  FontSizeLeg, 'FontSizeTicks',  FontSizeTicks, 'LegPos', 'bottom', 'Linespec', {'k-' , 'r-' , 'b-' , 'g-' , 'c-' , 'm-' },...
 'LineWidth',  2, 'Position', [ 246 28.65 106.875 31.9], 'SEM', 'on', 'Style', 'Classic', 'Tag', 'ERP_figure', 'Transparency',  0.8, 'xscale',...
 xscale, 'YDir', YDir, 'yscale', yscale );


%
% plot CTs
%
ERP = pop_ploterps( ERP,  BinArray_CT ,  ChanArray , 'Axsize', [ 0.05 0.08], 'BinNum', 'on', 'Blc', baseline, 'Box', [ 5 1], 'ChLabel', 'on',...
 'FontSizeChan',  FontSizeChan, 'FontSizeLeg',  FontSizeLeg, 'FontSizeTicks',  FontSizeTicks, 'LegPos', 'bottom', 'Linespec', {'k-' , 'r-' , 'b-' , 'g-' , 'c-' , 'm-' },...
 'LineWidth',  2, 'Position', [ 246 28.65 106.875 31.9], 'SEM', 'on', 'Style', 'Classic', 'Tag', 'ERP_figure', 'Transparency',  0.8, 'xscale',...
 xscale, 'YDir', YDir, 'yscale', yscale );

%
% plot difference: ct02-it01
%
ERP = pop_ploterps( ERP,  BinArray_DIFF ,  ChanArray , 'Axsize', [ 0.05 0.08], 'BinNum', 'on', 'Blc', baseline, 'Box', [ 5 1], 'ChLabel', 'on',...
 'FontSizeChan',  FontSizeChan, 'FontSizeLeg',  FontSizeLeg, 'FontSizeTicks',  FontSizeTicks, 'LegPos', 'bottom', 'Linespec', {'k-' , 'r-' , 'b-' , 'g-' , 'c-' , 'm-' },...
 'LineWidth',  2, 'Position', [ 246 28.65 106.875 31.9], 'SEM', 'on', 'Style', 'Classic', 'Tag', 'ERP_figure', 'Transparency',  0.8, 'xscale',...
 xscale, 'YDir', YDir, 'yscale', yscale );