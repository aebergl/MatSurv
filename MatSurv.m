function [varargout] = MatSurv(TimeVar, EventVar, GroupVar, varargin)
%MatSurv Survival analysis with Kaplan-Meier plot, log rank test and risk table
%
%
% USAGE:
% MatSurv(TimeVar, EventVar, GroupVar,'param', value, ...) creates a Kaplan-Meier plot
%   [p] = MatSurv( ... ) returns the log rank p-value
%   [p, fh] = MatSurv( ... ) returns both p-value and figure handle
%   [p, fh, stats] = MatSurv( ... ) returns additions stats from log rank test
%
% INPUTS:
% * 'TimeVar' is a vector with numeric time to event values of events, either observed or
%   censored. Values equal or less than zero will be removed by default
%
% * 'EventVar' is a vector or cell array defining events or censored
%   observation. Events are defined with a 1 and censored point with a 0. By
%   default 'Dead', 'Deceased', 'Relapsed', 'Yes' are considered as events. 
%   'Alive', 'Living', 'Not Relapsed', 'DiseaseFree', 'No' are considers as censored
%   'EventDefinition' can be used to define other types of events
%
% * 'GroupVar' is a vector or cell array that defines the differnt groups.
%   If it is continues variable median cut will be used as defining two
%   groups by default. 
%
% OUTPUTS:
% * p       : log rank p-value
% * fh      : figure handle to KM-plot figure
% * stats   : Additional statistics from the log rank test
%

p.addParameter('EventDefinition',[]);
p.addParameter('FlipGroupOrder',0);
p.addParameter('FlipColorOrder',0);
p.addParameter('BaseFontSize',16);

% OTHER PARAMETERS (passed as parameter-value pairs)
% * 'NoPlot': A true/false value which, if true, no figure is created
%   (defaut: false)

% * 'CalcP': A true/false value which, if true, a log rank test is
%   performed and diplayd on the KM-plot. (defaut: true)

% * 'CutPoint': Either a string or scalar/vector with cut points to be used
%   for defining groups baased on a continious 'GroupVar' input variable
%   Allowed names are: 'Median' or 'Quartile'
%   If a scalar or vector is used the groups will be defined based on the
%   cut points. (defaut: 'median')
%
% * 'GroupsToUse': Cell array defining what groups to use from the GroupVar
%   variable. Works only if GrouVar is a cell array. (default: all groups are used)
%
% * 'GroupOrder': A cell array defning the group order to be used in the
%   legend. (default: Groups are sorted alphabetically)
%
% * 'EventDefinition': Two element cell array where the first cell defines
%   the event and the second censored values. Example {'Dead,'Alive'}
%   (default = true)
%
% * 'NaNColor': A 3-element [R G B] vector specifying the color used to display NaN
%   or missing value. [0 0 0] corresponds to black and [1 1 1] to white. By
%   default MATLAB displays NaN values using the color assigned to the
%   lowest value in the colormap. Specifying this option automatically sets
%   the 'UseFigureColormap' option to false because the color mapping must
%   be computed prior to setting the nan color.
%
% * 'MinColorValue': A scalar number corresponding to the value of the data
%   that is mapped to the lowest color of the colormap. By default this is
%   the minimum value of the matrix input.
%
% * 'MaxColorValue': A scalar number corresponding to the value of the data
%   that is mapped to the highest color of the colormap. By default this is
%   the maximum value of the matrix input.
%
% * 'Parent': Handle to an axes object
%
% * 'TextColor': Either a color specification of all the text displayed on
%   the image or a string 'xor' which sets the EraseMode property of the text
%   objects to 'xor'. This will display all the text labels in a color that
%   contrasts its background.
%
% * 'FontSize': The initial fontSize of the text labels on the image. As
%   the image size is scaled the fontSize is shrunk appropriately.
%
% * 'ColorBar': Display colorbar. The corresponding value parameter should
%   be either logical 1 or 0 or a cell array of any additional parameters
%   you wish to pass to the colorbar function (such as location)
%
% * 'GridLines': Draw grid lines separating adjacent sections of the
%   heatmap. The value of the parameter is a LineStyle specification, for example,
%   :, -, -. or --. By default, no grid lines are drawn.
%
% * 'TickAngle': Angle of rotation of tick labels on x-axis. (Default: 0)
%
% * 'ShowAllTicks': Set to 1 or true to force all ticks and labels to be
%   drawn. This can make the axes labels look crowded. (Default: false)
%
% * 'TickFontSize': Font size of the X and Y tick labels. Default value is
%   the default axes font size, usually 10. Set to a lower value if many
%   tick labels are being displayed
%
% * 'TickTexInterpreter': Set to 1 or true to render tick labels using a TEX
%   interpreter. For example, '_b' and '^o' would be rendered as subscript
%   b and the degree symbol with the TEX interpreter. This parameter is only
%   available in MATLAB R2014b and above (Default: false)
%
% OUTPUTS:
% * hImage: handle to the image object
% * hText : handle to the text objects (empty if no text labels are drawn)
% * hTick : handle to the X-tick label text objects if tick angle is not 0
%           (empty otherwise)
%
% Notes:
% * The 'money' colormap displays a colormap where 0 values are mapped to
%   white, negative values displayed in varying shades of red and positive
%   values in varying shades of green
% * The 'red' colormap maps 0 values to white and higher values to red
%
% EXAMPLES:
% data = reshape(sort(randi(100,10)),10,10)-50;
% heatmap(data, cellstr(('A':'J')'), mean(data,2), '%0.0f%%',...
%         'Colormap', 'money', 'Colorbar', true, 'GridLines', ':',...
%         'TextColor', 'b')
% For detailed examples, see the associated document heatmap_examples.m

% Copyright The MathWorks, Inc. 2009-2014

% *** Anders Berglund ***


% Check TimeVar, EventVar, GroupVar variables
if nargin < 3
    error('MatSurv requires at least 3 input argument');
end

% Load test data
if isempty(TimeVar) && isempty(EventVar) && isempty(GroupVar)
    [TimeVar, EventVar, GroupVar] = SurvMatLoadTestData;
    varargin =[varargin,{'TimeUnit','Weeks'}];
end

% Check that they are all vectors
if min(size(TimeVar)) ~= 1 || min(size(EventVar)) ~= 1 || min(size(GroupVar)) ~= 1
    error('TimeVar, EventVar, GroupVar must all be vectors or cell arrays');
end

% Check that they all are equal length
if (length(TimeVar) ~= length(EventVar)) || (length(TimeVar) ~= length(GroupVar)) || (length(EventVar) ~= length(GroupVar))
    error('TimeVar, EventVar, GroupVar must have equal length');
end

%Parse input and set default values
options = MatSurvParseInput(varargin{:});

% Check input and clean input data
[TimeVar, EventVar, GroupVar] = MatSurvCleanData(TimeVar, EventVar, GroupVar, options);

% Define events 1=event, 0=no event but it also checks for dead/alive etc use
% EventDefinition parameter for full control
[EventVarBin] = MatSurvDefineEventVar(EventVar, options);

% Censor data if TimeMax is given
if ~isempty(options.TimeMax)
    [TimeVar, EventVarBin] = MatSurvCensorTimeMax(TimeVar, EventVarBin, options);
end

% CreatGroups based on GroupVar and create DATA structure
[DATA] = MatSurvCreateGroups(TimeVar, EventVarBin, GroupVar, options);

% Flip Group Ordering
if options.FlipGroupOrder
    DATA.GROUPS = DATA.GROUPS(DATA.numGroups:-1:1);
end

% Creat Survival table for plotting
[DATA] = MatSurvCreateTable(DATA);

% Do log rank test
if options.CalcP
    [p,stats] = MatSurvLogRank(DATA);
else
    p=[];
    stats=[];
end

if ~options.NoPlot
    % Create Figure Window
    fh=figure('Name','MatSurv KM-Plot','Color','w','Tag','MatSurv KM-Plot figure');
    
    %Create Axes
    axh_KM = axes(fh,'Position',options.KM_position,'NextPlot','add','tag','KM-Plot');
    axh_RT = axes(fh,'Position',options.RT_position,'tag','Risk Table');
    % No axis for the Risk Table
    axh_RT.XAxis.Visible='off';
    axh_RT.YAxis.Visible='off';
    
    % Adjust Colors for user input
    if ischar(options.LineColor)
        cMAP = colormap(options.LineColor);
    elseif ismatrix(options.LineColor)
        cMAP = options.LineStyle;
    end
    cMAP = cMAP(1:DATA.numGroups,:);
    if options.FlipColorOrder
        cMAP = flipud(cMAP);
    end
    % Adjust line style
    if ischar(options.LineStyle)
        LineStyles = cell(DATA.numGroups,1);
        LineStyles(:) = {options.LineStyle};
    elseif numel(options.LineStyle) == 1
        LineStyles = cell(DATA.numGroups,1);
        LineStyles(:) = options.LineStyle;
    elseif iscell(options.LineStyle)
        LineStyles = options.LineStyle;
    end
    % Adjust censoring markers
    if ischar(options.CensorLineColor) && strcmpi('same',options.CensorLineColor)
        cMAPCensor = cMAP;
    elseif ismatrix(options.CensorLineColor)
        cMAPCensor = options.CensorLineColor;
    end
    
    % Create stairs
    S=ones(DATA.numGroups,1);
    for i=1:DATA.numGroups
        S(i)=stairs(axh_KM,DATA.GROUPS(i).KM_ALL(:,1),DATA.GROUPS(i).KM_ALL(:,2),'Color',cMAP(i,:),'Linewidth',options.LineWidth,'LineStyle',LineStyles{i});
        if ~isempty(DATA.GROUPS(i).Censored_Points)
            % Draw marks for censored points
            line(axh_KM,[DATA.GROUPS(i).Censored_Points(:,1)'; DATA.GROUPS(i).Censored_Points(:,1)'],...
                [DATA.GROUPS(i).Censored_Points(:,2)'-options.CensorLineLength ; DATA.GROUPS(i).Censored_Points(:,2)'+options.CensorLineLength],...
                'Color',cMAPCensor(i,:),'Linewidth',options.CensorLineWidth);
        end
    end
    
    %Fix Y-Axis
    % Limir range from 0 to 1
    axh_KM.YLim = [0 1];
    axh_KM.YTick = options.YTick;
    YMinorStep =  (options.YTick(2) - options.YTick(1) ) / (1+options.YMinorTick);
    axh_KM.YAxis.MinorTickValues = YMinorStep:YMinorStep:1;
    axh_KM.YAxis.MinorTick = 'on';
    axh_KM.YAxis.TickDirection = 'out';
    
    % Y label
    axh_KM.YAxis.FontSize=options.YTickFontSize;
    ylabel(axh_KM,options.Ylabel,'FontSize',options.YLabelFontSize,options.YlabelOptions{:});
    
    % X label
    axh_KM.XAxis.FontSize=options.XTickFontSize;
    if isempty(options.Xlabel)
        xlabel_str = sprintf('Time (%s)',options.TimeUnit);
    else
        xlabel_str =options.Xlabel;
    end
    xlabel(axh_KM,xlabel_str,'FontSize',options.XLabelFontSize,options.XlabelOptions{:});
    axh_KM.XAxis.TickDirection = 'out';
    
    % Title
    if ~isempty(options.Title)
        title(axh_KM,options.Title,'FontSize',18,options.TitleOptions{:});
    end
    
    % Set legend
    h_LE=legend(S,[DATA.GROUPS(:).GroupName]);
    h_LE.Box='off';
    h_LE.FontSize=options.LegendFontSize;
    
    if ~isempty(options.XLim)
        axh_KM.XLim = [0 options.XLim];
    end
    
    % Get Xticks
    %max_X = max(max(DATA.GROUPS.TimeVar));
    max_X = axh_KM.XLim(2);
    Nudge_X = max_X / 50;
    
    if ~isempty(options.Xstep)
        axh_KM.XTick = 0:options.Xstep:max_X;
    end
    if ~isempty(options.XTicks)
        axh_KM.XTick = options.XTicks;
    end
    axh_RT.XTick=axh_KM.XTick;
    axh_KM.XAxis.MinorTick = 'on';
    XMinorStep =  (axh_KM.XTick(2) - axh_KM.XTick(1) ) / (1+options.XMinorTick);
    axh_KM.XAxis.MinorTickValues = XMinorStep:XMinorStep:axh_KM.XTick(end);
    axh_KM.LineWidth = 1.5;
    
    if options.CalcP
        txt_str = sprintf('p = %.3g',p);
        text(axh_KM,Nudge_X,0.1,txt_str,'FontSize',options.PvalFontSize,'tag','p-value')
    end
    
    % And now to the Risk table
    
    % Get number of samples for each time point
    RT_X = zeros(length(axh_KM.XTick),DATA.numGroups);
    for i = 1:length(axh_KM.XTick)
        for j = 1:DATA.numGroups
            %RT_X(i,j) = sum(DATA.GROUPS(j).TimeVar > axh_KM.XTick(i) & DATA.GROUPS(j).EventVar == 1) + sum(DATA.GROUPS(j).TimeVar >= axh_KM.XTick(i) & DATA.GROUPS(j).EventVar == 0);
            RT_X(i,j) = sum(DATA.GROUPS(j).TimeVar >= axh_KM.XTick(i));
        end
        
    end
    axh_RT.YLim = [0.5 DATA.numGroups + 0.5];
    axh_RT.YTick = 1:DATA.numGroups;
    linkaxes([axh_RT,axh_KM],'x')
    
    % Color OptionsFor Risk Table
    if ischar(options.RT_Color) && strcmpi('same',options.RT_Color)
        cMAP_RT = cMAP;
    elseif ismatrix(options.RT_Color)
        cMAP_RT = options.RT_Color;
        cMAP_RT = repmat(cMAP_RT,DATA.numGroups,1);
    end
    
    for i = 1:length(axh_KM.XTick)
        for j = 1:DATA.numGroups
            %sprintf('%u',RT_X(i,j))
            text(axh_RT,axh_RT.XTick(i),axh_RT.YTick(end-j+1),sprintf('%u',RT_X(i,j)),...
                'HorizontalAlignment','center','VerticalAlignment','middle',...
                'FontSize',options.RT_FontSize,'Color',cMAP_RT(j,:))
        end
    end
    if options.RT_YLabel
        left_pos = axh_RT.Children(end).Extent(1);
        nudge_x = abs(axh_RT.XLim(2) - axh_RT.XLim(1))/100;
        
        line(axh_RT,[left_pos-nudge_x left_pos-nudge_x],[axh_RT.YTick(1)-0.5 axh_RT.YTick(end)+0.5],'color','k','clipping','off','LineWidth',1.25)
        for j = 1:DATA.numGroups
            text(axh_RT,left_pos-(nudge_x*2),axh_RT.YTick(end-j+1),DATA.GROUPS(j).GroupName,...
                'HorizontalAlignment','right','VerticalAlignment','middle',...
                'FontSize',options.RT_FontSize,'Color',cMAP_RT(j,:),'FontWeight','bold')
        end
    end
    % Title
    if ~isempty(options.RT_Title)
        ht = title(axh_RT,options.RT_Title,'FontSize',14,options.TitleOptions{:});
        ht.VerticalAlignment='middle';
    end
else
    fh = [];
end

if nargout > 0
    varargout{1} = p;
end
if nargout > 1
    varargout{2} = fh;
end
if nargout > 2
    varargout{3} = stats;
end


end

function params = MatSurvParseInput(varargin)
%Parse input and set defualt values
p = inputParser;
p.addParameter('NoPlot',false);
p.addParameter('CalcP',1);
p.addParameter('CutPoint','Median');
p.addParameter('GroupOrder',[]);
p.addParameter('GroupsToUse',[]);
p.addParameter('EventDefinition',[]);
p.addParameter('FlipGroupOrder',0);
p.addParameter('FlipColorOrder',0);
p.addParameter('BaseFontSize',16);

p.addParameter('NoWarnings',false);
p.addParameter('TimeMin',0, @(x)isnumeric(x) && isscalar(x));
p.addParameter('TimeMax',[], @(x)isnumeric(x) && isscalar(x));
p.addParameter('TimeUnit','Months');

% Figure Options
p.addParameter('KM_position',[0.2 0.4 0.75 0.45]);
p.addParameter('RT_position',[0.2 0.05 0.75 0.20]);

% KM plot options
p.addParameter('Xstep',[], @(x)isnumeric(x) && isscalar(x));
p.addParameter('XTicks',[], @(x)isnumeric(x) && isvector(x));
p.addParameter('XMinorTick',1, @(x)isnumeric(x) && isscalar(x));

p.addParameter('XLim',[], @(x)isnumeric(x) && isscalar(x));
p.addParameter('LineColor','lines');
p.addParameter('LineWidth',2);
p.addParameter('LineStyle','-');
p.addParameter('CensorLineWidth',2);
p.addParameter('CensorLineLength',0.02);
p.addParameter('CensorLineColor','same');

p.addParameter('Xlabel',[]);
p.addParameter('XlabelOptions',cell(0,0));
p.addParameter('XLabelFontSize',16);
p.addParameter('XTickFontSize',14);

p.addParameter('Ylabel','Survival Probability');
p.addParameter('YlabelOptions',cell(0,0));
p.addParameter('YLabelFontSize',16);
p.addParameter('YTickFontSize',14);
p.addParameter('YTick',0:0.2:1);
p.addParameter('YMinorTick',1);

p.addParameter('Title',[]);
p.addParameter('TitleOptions',cell(0,0));
p.addParameter('LegendFontSize',14);
p.addParameter('PvalFontSize',16);

% Risk table plot options
p.addParameter('RT_FontSize',16);
p.addParameter('RT_Color','same');
p.addParameter('RT_YLabel',1);
p.addParameter('RT_Title',[]);
p.addParameter('RT_TitleOptions',cell(0,0));

parse(p,varargin{:});
params = p.Results;

end

function [p,stats] = MatSurvLogRank(DATA)

% Merge tables from all groups
KM_ALL = vertcat(DATA.GROUPS.KM_Events);

% Get all time points with events
tf = KM_ALL(:,1);
tf = unique(tf);

% allocate matrices
n = length(tf);
mf = zeros(n,DATA.numGroups);
nf = zeros(n,DATA.numGroups);
ef = zeros(n,DATA.numGroups-1);

% Assign values
for i = 1:DATA.numGroups
    % Need to add time entries including censored ones
    tf_in = unique([tf;DATA.GROUPS(i).TimeVar]);
    [KM_Events, ~, ~] = MatSurvCalculateTables(tf_in,DATA.GROUPS(i).TimeVar,DATA.GROUPS(i).EventVar,tf);
    nf(:,i) = KM_Events(:,2);
    mf(:,i) = KM_Events(:,3);
    
end

nf_sum = sum(nf,2);
mf_sum = sum(mf,2);
for i = 1:DATA.numGroups-1
    ef(:,i) = (nf(:,i)  ./ nf_sum) .* mf_sum;
end
d = sum(mf(:,1:end-1)-ef)';
Var_OE=zeros(n,DATA.numGroups-1);
for i = 1:DATA.numGroups-1
    Var_OE(:,i) = (nf(:,i) .* (nf_sum - nf(:,i)) .* mf_sum .*(nf_sum - mf_sum)) ./ (nf_sum.^2 .* (nf_sum -1));
    %Var_OE(:,i) = (nf(:,i) .* (nf_sum - nf(:,i)) .* mf(:,i) .*(nf_sum - mf_sum)) ./ (nf_sum.^2 .* (nf_sum -1));
end
Var_OE_sum =sum(Var_OE);
Cov_OE = zeros(n,(DATA.numGroups-1)*(DATA.numGroups-2)/2);
if DATA.numGroups > 2
    counter = 0;
    for i = 1:DATA.numGroups-2
        for j = i+1:DATA.numGroups-1
            counter = counter + 1;
            Cov_OE(:,counter) = ( -nf(:,i) .* nf(:,j) .* mf_sum .* (nf_sum - mf_sum)) ./ (nf_sum.^2 .* (nf_sum -1));
        end
    end

    Cov_OE_sum = sum(Cov_OE);
    V = zeros(DATA.numGroups-1);
    V(tril(true(DATA.numGroups-1),-1))=Cov_OE_sum;
    V(~tril(true(DATA.numGroups-1),0))=Cov_OE_sum;
    V(1:size(V,1)+1:end) = Var_OE_sum;
else
    V = Var_OE_sum;
end
LogRank_Stat = d'/V*d;
p = 1 - gammainc(LogRank_Stat/2,(DATA.numGroups-1)/2);
stats.Chi2 = LogRank_Stat;

% elseif DATA.numGroups ==2
%     ef = zeros(n,DATA.numGroups);
%
%     npn = nf(:,1) + nf(:,2);
%     mpm = mf(:,1) + mf(:,2);
%     ntn = nf(:,1) .* nf(:,2);
%     for i = 1:DATA.numGroups
%         ef(:,i) = (nf(:,i)  ./ npn) .* mpm;
%     end
%     O_E = mf-ef;
%     sum_OE = sum(O_E)
%     Var_O_E = (ntn .* mpm .* (npn - mpm)) ./ ( npn.^2 .* (npn - 1));
%     sum_Var_O_E =sum(Var_O_E)
%     LogRank_Stat = sum_OE(1)^2/sum_Var_O_E
%     p = 1-chi2cdf(LogRank_Stat,1);
% end
end


function [DATA] = MatSurvCreateTable(DATA)

for i=1:DATA.numGroups
    % Get unique time points including censored and add a leading 0
    tf = [0; unique(DATA.GROUPS(i).TimeVar)];
    
    [KM_Events, KM_ALL, Censored_Points] = MatSurvCalculateTables(tf,DATA.GROUPS(i).TimeVar,DATA.GROUPS(i).EventVar,[]);
    
    DATA.GROUPS(i).KM_Events=KM_Events;
    DATA.GROUPS(i).KM_ALL=KM_ALL;
    DATA.GROUPS(i).Censored_Points=Censored_Points;
end

end

function [KM_Events, KM_ALL, Censored_Points] = MatSurvCalculateTables(tf,TimeVar,EventVar,tf_out)

% Calculate number  of samples for each time point including censored
mf = sum(repmat(TimeVar,1,length(tf)) == repmat(tf',length(TimeVar),1))';

%Calculate number of samples left at each time point
mf_cumsum = cumsum(mf);
nf = ones(length(tf),1) * length(TimeVar);
nf(2:end) = nf(2:end) - mf_cumsum(1:end-1);

% Find censored points
indx_censor = (EventVar == 0);
tfq = unique(TimeVar(indx_censor));
mfq = sum(repmat(TimeVar(indx_censor),1,length(tfq)) == repmat(tfq',length(TimeVar(indx_censor)),1))' ;

% Find time points where there are censored data
[~,tf_indx,~]=intersect(tf,tfq,'stable');

%Adjust counts for censored data, a 0 will mean that only censored data
%was observed at that time
mf_true = mf;
mf_true(tf_indx) = mf((tf_indx)) - mfq;

% Calculate fraction alive
S = cumprod(1-(mf_true./nf));

% get index of censored samples for plotting
indx_censored = ((mf_true - mf) < 0);

%
if isempty(tf_out) % function called for plotting
    indx_observed = (mf_true ~= 0);
else % fynction called for log rank test
    [~,~,indx_observed] = intersect(tf_out,tf);
end

KM_Events = [tf(indx_observed) nf(indx_observed) mf_true(indx_observed)];
KM_ALL = [tf S nf];
Censored_Points=[tf(indx_censored) S(indx_censored)];

end
function [DATA] = MatSurvCreateGroups(TimeVar, EventVarBin, GroupVar, options)
% Create Group structure
DATA.numGroups = 0;
DATA.GROUPS = struct('GroupName',{},'TimeVar',[],'EventVar',[]);

% Define set of Groups to use
if ~isempty(options.GroupsToUse) % User defined Groups to use
    DATA.numGroups = numel(options.GroupsToUse);
    for i = 1:DATA.numGroups
        if iscell(GroupVar)
            indx_group = strcmp(options.GroupsToUse(i),GroupVar);
            DATA.GROUPS(i).GroupName = options.GroupsToUse(i);
        elseif isnumeric(GroupVar)
            indx_group = (options.GroupsToUse(i) == GroupVar);
            DATA.GROUPS(i).GroupName = {num2str(options.GroupsToUse(i))};
        end
        DATA.GROUPS(i).TimeVar = TimeVar(indx_group);
        DATA.GROUPS(i).EventVar = EventVarBin(indx_group);
    end
    
    % If the Groupvariable is a cell vector
elseif iscell(GroupVar)
    Unique_Groups = unique(GroupVar);
    DATA.numGroups = length(Unique_Groups);
    for i = 1:DATA.numGroups
        indx_group = strcmp(Unique_Groups(i),GroupVar);
        DATA.GROUPS(i).GroupName = Unique_Groups(i);
        DATA.GROUPS(i).TimeVar = TimeVar(indx_group);
        DATA.GROUPS(i).EventVar = EventVarBin(indx_group);
    end
    % If the Groupvariable is a numerical vector
elseif (strcmpi('Median',options.CutPoint) || isscalar(options.CutPoint)) && isnumeric(GroupVar)
    if strcmpi('Median',options.CutPoint)
        Cut_Val = median(GroupVar);
        tmp_txt = 'Median';
    elseif isscalar(options.CutPoint)
        Cut_Val = options.CutPoint;
        tmp_txt = 'Fixed value';
    end
    DATA.numGroups = 2;
    indx_Above  = (GroupVar > Cut_Val);
    indx_Below  = (GroupVar <= Cut_Val);
    
    DATA.GROUPS(1).GroupName = {sprintf('%s > %g',tmp_txt,Cut_Val)};
    DATA.GROUPS(1).TimeVar = TimeVar(indx_Above);
    DATA.GROUPS(1).EventVar = EventVarBin(indx_Above);
    
    DATA.GROUPS(2).GroupName = {sprintf('%s <= %g',tmp_txt,Cut_Val)};
    DATA.GROUPS(2).TimeVar = TimeVar(indx_Below);
    DATA.GROUPS(2).EventVar = EventVarBin(indx_Below);
    
elseif strcmpi('Quartile',options.CutPoint)  && isnumeric(GroupVar)
    Cut_Val = prctile(GroupVar,[25 75]);
    tmp_txt = 'Quartile';
    indx_Above = (GroupVar >= Cut_Val);
    indx_Below  = (GroupVar <= Cut_Val);
    DATA.numGroups = 2;
    DATA.GROUPS(1).GroupName = {sprintf('%s > %g',tmp_txt,Cut_Val)};
    DATA.GROUPS(1).TimeVar = TimeVar(indx_Above);
    DATA.GROUPS(1).EventVar = EventVarBin(indx_Above);
    DATA.GROUPS(2).GroupName = {sprintf('%s <= %g',tmp_txt,Cut_Val)};
    DATA.GROUPS(2).TimeVar = TimeVar(indx_Below);
    DATA.GROUPS(2).EventVar = EventVarBin(indx_Below);
end

end

function [TimeVar, EventVarBin] = MatSurvCensorTimeMax(TimeVar, EventVarBin, options)
indx_TimeMax = (TimeVar > options.TimeMax);
TimeVar(indx_TimeMax) = TimeMax;
EventVarBin(indx_TimeMax) = 0;

end

function [EventVarBin] = MatSurvDefineEventVar(EventVar, options)

% Set all entries to zeros
EventVarBin = zeros(size(EventVar));

if islogical(EventVar) % Set TRUE to 1
    EventVarBin(EventVar) = 1;
elseif isnumeric(EventVar) % set ones to 1
    EventVarBin(EventVar == 1) = 1;
elseif iscell(EventVar)
    if ~isempty(options.EventDefinition) % Set values based on user input
        indx_Event = strcmp(options.EventDefinition{1},EventVar);
        indx_NoEvent = strcmp(options.EventDefinition{2},EventVar);
        if sum(indx_Event) + sum(indx_NoEvent) == length(EventVar)
            EventVarBin(indx_Event) = 1;
            EventVarBin(indx_NoEvent) = 0;
        else
            error('Event variable do not match event type defined in options.EventDefinition')
        end
    else % Set values based on common event types such as dead/alive
        indx_Event = strcmpi('dead',EventVar) | strcmpi('DECEASED',EventVar) | strcmpi('Relapsed',EventVar) |  strcmpi('Yes',EventVar) ;
        indx_NoEvent = strcmpi('alive',EventVar) | strcmpi('LIVING',EventVar) | strcmpi('NotRelapsed',EventVar) | strcmpi('DiseaseFree',EventVar)| strcmpi('No',EventVar);
        if sum(indx_Event) + sum(indx_NoEvent) == length(EventVar)
            EventVarBin(indx_Event) = 1;
            EventVarBin(indx_NoEvent) = 0;
        else
            error('Event variable has non recognazed type. please check EventVar')
        end
    end
else
    error('Non supported Event variable input')
end

end

function [TimeVar, EventVar, GroupVar] = MatSurvCleanData(TimeVar, EventVar, GroupVar, options)
% Functions to check and cleanup inout data
%
% Make sure that TimeVar, EventVar, GroupVar are column vectors and not row
% vectors
if size(TimeVar,1) == 1
    TimeVar = TimeVar';
end
if size(EventVar,1) == 1
    EventVar = EventVar';
end
if size(GroupVar,1) == 1
    GroupVar = GroupVar';
end


% Check time variable for missing data and timepoints < TimeMin
rem_indx_time = ( isnan(TimeVar) | (TimeVar <= options.TimeMin) );

% Check time variable for missing data and/or empty cells and cells with NA
if isnumeric(EventVar) || islogical(EventVar)
    rem_indx_event = isnan(EventVar);
elseif iscell(EventVar)
    rem_indx_event = any(isempty(EventVar), strcmpi('NA',EventVar));
end

% Check time variable for missing data and/or empty cells and cells with NA
if isnumeric(GroupVar)
    rem_indx_group = isnan(GroupVar);
elseif iscell(GroupVar)
    rem_indx_group = (isempty(GroupVar) | strcmpi('NA',GroupVar));
end

% Merge all indexes

rem_indx = (rem_indx_time | rem_indx_event |rem_indx_group);

if sum(rem_indx) > 0
    TimeVar(rem_indx) = [];
    EventVar(rem_indx) = [];
    GroupVar(rem_indx) = [];
    if ~options.NoWarnings
        fprintf('*********************************\n');
        fprintf('\n');
        fprintf('%u sample have been removed',sum(rem_indx));
        fprintf('removed samples had missind data or time < %g\n',options.TimeMin)
    end
    
end

numEvenTypes = unique(EventVar);

if numEvenTypes > 2
    error('More then 2 event types in the Event variable');
end

end



function [TimeVar, EventVar, GroupVar] = SurvMatLoadTestData

t1=[6 6 6 7 10 13 16 22 23 6 9 10 11 17 19 20 25 32 32 34 35]';
t2=[1 1 2 2 3 4 4 5 5 8 8 8 8 11 11 12 12 15 17 22 23]';
TimeVar=[t1;t2];

e1=[1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0]';
e2=ones(21,1);
EventVar = [e1;e2];

g1=cell(size(t1));
g1(:) = {'Group 1'};

g2=cell(size(t2));
g2(:) = {'Group 2'};

GroupVar = [g1;g2];

end
