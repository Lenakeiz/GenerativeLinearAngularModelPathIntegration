%% Script to create output for Fig. 1 - parameter comparisons between Young, healthy Elderly, MCI merged 
% Andrea Castegnaro, UCL, 2022 uceeaca@ucl.ac.uk
% Fits the model on all of the groups and run a two-way Anova
% (group*condition) on Young Controls, Elderly Controls and MCI negative
% Output: for each parameter fitted by the model output one boxplot with
% the three groups and performance splitted by environmental condition and
% a boxplot with three groups performance averaged across environmental
% conditions

% Preparing the data
GLAMPI_PrepareBaseConfig;

% Preprocessing the data
GLAMPI_PreprocessData;

% Model fitting
config.ModelName        =   "beta_k_g2_g3_m3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "m3", "sigma", "nu"];

config.NumParams        =   length(config.ParamName);

GLAMPI;

% Generating color scheme for our paper
ColorPattern;

%% Preparing output
config.ResultFolder = pwd + "/Output/Fig1/MCI_Only";
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Collecting information from output
[MCINegPropDist, MCINegPropAng]                   = getProportionalLinearAndAngularError(MCINeg, MCINeg.Results.estimatedParams);
[MCIPosPropDist, MCIPosPropAng]                   = getProportionalLinearAndAngularError(MCIPos, MCIPos.Results.estimatedParams);

% TwowayAnova Analysis
config.type = "ProportionalDistance";
[anova_tab_dist,multicomp_tab1_dist,multicomp_tab2_dist, multicomp_tab12_dist]     = TwowayAnova_Behavioural_MCIOnly(MCIPosPropDist, MCINegPropDist, config);

% TwowayAnova Analysis
config.type = "ProportionalAngle";
[anova_tab_angle,multicomp_tab1_angle,multicomp_tab2_angle, multicomp_tab12_angle] = TwowayAnova_Behavioural_MCIOnly(MCIPosPropAng, MCINegPropAng, config);

% Display means and std for each group
disp("%%%%%%%%% Proportional Distance Error %%%%%%%%%");
disp(["MCI Positive : " num2str(mean(mean(MCIPosPropDist,1,"omitnan"),"omitnan")) " +- " std(mean(MCIPosPropDist,1,"omitnan"),"omitnan")]);
disp(["MCI Negative : " num2str(mean(mean(MCINegPropDist,1,"omitnan"),"omitnan")) " +- " std(mean(MCINegPropDist,1,"omitnan"),"omitnan")]);

TTest_NominalValue(mean(MCIPosPropDist,1,"omitnan")',1,'MCI positive', 'both');
TTest_NominalValue(mean(MCINegPropDist,1,"omitnan")',1,'MCI negative', 'both');

disp("%%%%%%%%% Proportional Angular Error %%%%%%%%%");
disp(["MCI Positive : " num2str(mean(mean(MCIPosPropAng,1,"omitnan"),"omitnan")) " +- " std(mean(MCIPosPropAng,1,"omitnan"),"omitnan")]);
disp(["MCI Negative : " num2str(mean(mean(MCINegPropAng,1,"omitnan"),"omitnan")) " +- " std(mean(MCINegPropAng,1,"omitnan"),"omitnan")]);

TTest_NominalValue(mean(MCIPosPropAng, 1, "omitnan")', 1, 'MCI positive', 'both');
TTest_NominalValue(mean(MCINegPropAng, 1, "omitnan")', 1, 'MCI negative', 'both');

%% Proportional Distance Error
plotInfo.defaultTextSize = 20;
plotInfo.defaultLineSize = 2;
plotInfo.barFaceAlpha = 0.5;
plotInfo.barLineWidth = 0.3;
plotInfo.barWidth = 0.5;
plotInfo.scatterFaceAlpha = 0.5;
plotInfo.scatterEdgeAlpha = 0.8;
plotInfo.scatterDataSize = 32;
plotInfo.errorBarWidth = 2.0;
plotInfo.visible = "on";
plotInfo.dimensions = [100 100 250 250];

plotInfo.type = "ProportionalDistance";
plotInfo.YLabel = "Actual distance / correct distance";
plotInfo.yLim = [0 1.75];
plotInfo.ticksStep = 0.5;

BoxPlotOfFittedParamMergeCondition(MCIPosPropDist, MCINegPropDist, anova_tab_dist, multicomp_tab1_dist, config, plotInfo);

% Proportional Angular Error
plotInfo.defaultTextSize = 20;
plotInfo.defaultLineSize = 2;
plotInfo.barFaceAlpha = 0.5;
plotInfo.barLineWidth = 0.3;
plotInfo.scatterFaceAlpha = 0.5;
plotInfo.scatterEdgeAlpha = 0.8;
plotInfo.scatterDataSize = 32;
plotInfo.errorBarWidth = 2.0;
plotInfo.type = "ProportionalAngle";
plotInfo.YLabel = "Actual angle / correct angle";

plotInfo.yLim = [0 2.5];
BoxPlotOfFittedParamMergeCondition(MCIPosPropAng, MCINegPropAng, anova_tab_angle, multicomp_tab1_dist, config, plotInfo);

% Final cleanup to leave workspace as the end of the Preprocessing stage.
% Remove if you want to take a look at the output data.
% clearvars -except config YoungControls HealthyControls MCINeg MCIPos MCIUnk anova_tab_angle anova_tab_dist multicomp_tab1_dist multicomp_tab2_dist multicomp_tab12_dist multicomp_tab1_angle multicomp_tab2_angle multicomp_tab12_angle

%% ---------------------------------------------------------------------
function BoxPlotOfFittedParamMergeCondition(MCIPos, MCINeg, ANOVA_tab, ANOVA_multcomp_group, config, plotInfo)

MCIPosMean = mean(MCIPos, 1, "omitnan")';
MCINegMean = mean(MCINeg, 1, "omitnan")';

f = figure('visible', plotInfo.visible,'Position', plotInfo.dimensions);

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12)

colorForMCIPos = config.color_scheme_npg(6,:);
colorForMCINeg = config.color_scheme_npg(3,:);

% parameters set for controlling visual output
whisker_value               =   1.5;
box_lineWidth               =   plotInfo.barLineWidth;
box_widths_value            =   plotInfo.barWidth;
box_color_transparency      =   plotInfo.barFaceAlpha; %faceAlpha
median_lineWidth            =   2;
median_color                =   'k';
scatter_jitter_value        =   0.4;
scatter_markerSize          =   plotInfo.scatterDataSize;
scatter_marker_edgeColor    =   'k';
scatter_marker_edgeWidth    =   0.5;
scatter_color_transparency  =   plotInfo.scatterFaceAlpha; %faceAlpha
mean_scatter_multiplier     =   2.0;

hold on
%% Boxplot for each column in MCIPos
bp1 = boxplot(MCIPosMean, ...
    'Whisker',whisker_value, ...
    'symbol','', ... %symbol ='' making outlier invisible
    'Color','k', ...
    'Notch','on', ...
    'widths',box_widths_value,...
    'positions', 1);
set(bp1,'linewidth',box_lineWidth);

hold on
%% Boxplot for each column in MCI Neg
bp2 = boxplot(MCINegMean, ...
    'Whisker',whisker_value, ...
    'symbol','', ... %symbol ='' making outlier invisible
    'Color','k', ...
    'Notch','on', ...
    'widths',box_widths_value,...
    'positions', 2);
set(bp2,'linewidth',box_lineWidth);

%% Boxplot visual changes
% matlab has a lifo system for these
h = findobj(gca,'Tag','Box');
%get the MCI pos box
patch(get(h(2),'XData'),get(h(2),'YData'),colorForMCIPos,'FaceAlpha',box_color_transparency);
%get the MCI neg box
patch(get(h(1),'XData'),get(h(1),'YData'),colorForMCINeg,'FaceAlpha',box_color_transparency);

%% Median visual change
h=findobj(gca,'tag','Median');
for i = 1:length(h)
    h(i).LineWidth = median_lineWidth;
    h(i).Color = median_color;
end

%% Scatter plot for data and mean (MCIPos)
num_points = size(MCIPosMean,1);
hold on
x = 1*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
scatter(x, MCIPosMean, scatter_markerSize, ...
    'filled', ...
    'o', ...
    'MarkerEdgeColor',scatter_marker_edgeColor, ...
    'MarkerFaceColor',colorForMCIPos, ...
    'MarkerFaceAlpha',scatter_color_transparency,...
    'LineWidth',scatter_marker_edgeWidth);

hold on
%add errorbar
mean_Young = mean(MCIPosMean, "omitnan");
sem_Young = std(MCIPosMean,"omitnan")./sqrt(length(MCIPosMean(~isnan(MCIPosMean))));
errorbar(1,mean_Young,sem_Young,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18);
hold on
%add mean point
scatter(1, mean_Young, mean_scatter_multiplier*scatter_markerSize, 'd',...
    'filled','MarkerEdgeColor','k', ...
    'MarkerFaceColor','w', ...
    'LineWidth',scatter_marker_edgeWidth);

%% Scatter plot for data and mean (MCI neg)
num_points = length(MCINegMean);
hold on
x = 2*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
scatter(x, MCINegMean, scatter_markerSize, ...
    'filled', ...
    'o', ...
    'MarkerEdgeColor',scatter_marker_edgeColor, ...
    'MarkerFaceColor',colorForMCINeg, ...
    'MarkerFaceAlpha',scatter_color_transparency,...
    'LineWidth',scatter_marker_edgeWidth);

%add errorbar
mean_Hold = mean(MCINegMean, "omitnan");
sem_Hold = std(MCINegMean, "omitnan")./sqrt(length(MCINegMean(~isnan(MCINegMean))));
errorbar(2,mean_Hold,sem_Hold,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18);
hold on
%add mean point
scatter(2, mean_Hold, mean_scatter_multiplier*scatter_markerSize, 'd',...
    'filled','MarkerEdgeColor','k', ...
    'MarkerFaceColor','w', ...
    'LineWidth',scatter_marker_edgeWidth);

hold on
%add horizontal line
yline(1,Color=[1.0, 0.5, 0.0],LineStyle='--',LineWidth=2);

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'XColor'      , [.0 .0 .0], ...
    'YColor'      , [.0 .0 .0], ...
    'XTick'       , (1:2),...
    'XLim'        , [0.5, 2.5],...
    'XTickLabel'  , {'MCI+','MCI-'},...
    'LineWidth'   , 1.0        );

ylabel(plotInfo.YLabel);
ylim(plotInfo.yLim);
yticks(0:plotInfo.ticksStep:plotInfo.yLim(2));

% Extract pvalues from multiple comparison of group effect. Adding
% it to the figure
multicomp_result = ANOVA_multcomp_group;
% 1       2             3      
% Young   HealthyOld    MCI
PvalueMCIPosvsMCINeg = multicomp_result(1,6); % MCIPos vs. MCINeg

%% Add significance bars
AllP = [PvalueMCIPosvsMCINeg];
Xval = [[1,2]];
%select those P value smaller than 0.05 (only add line when p<0.05)
% * represents p<=0.05
% ** represents p<=1E-2
% *** represents p<=1E-3
PsigInd = AllP<0.05;
if sum(AllP<0.05)>0
    Xval_select = Xval(PsigInd,:);
    AllP_select = AllP(PsigInd);
    XXX = {};
    for i=1:sum(AllP<0.05)
        XXX{i} = Xval_select(i,:);
    end
    H=adjustablesigstar(XXX,AllP_select);
end

hold off

%% export figure
exportgraphics(f,config.ResultFolder+"/MergeCondsBox_"+plotInfo.type+".png",'Resolution',300);
exportgraphics(f,config.ResultFolder+"/MergeCondsBox_"+plotInfo.type+".pdf",'Resolution',300, 'ContentType','vector');

end

%% --------------------------------------------------------------------- 
function [PropDist, PropAng] = getProportionalLinearAndAngularError(Group, Glampi_data)
% Get the proportiona distance and angular error for all particpants.
% Each value is the mean value per participant per environmental condition

    PropDist = [];
    PropAng  = [];

    for TRIAL_FILTER = 1:3
        
        CondPropDist  = Group.Results.PropDistErr{TRIAL_FILTER};
        CondPropAng   = Group.Results.PropAngErr{TRIAL_FILTER};
        subjectSize   = size(CondPropDist,2);

        CondPropDistMeanSubjs  = zeros(1,subjectSize);
        CondPropAngleMeanSubjs = zeros(1,subjectSize);
        
        for subj=1:subjectSize
            CondPropDistMeanSubjs(1,subj)  = mean(cell2mat(CondPropDist{1,subj}),"omitnan");
            CondPropAngleMeanSubjs(1,subj) = mean(cell2mat(CondPropAng{1,subj}),"omitnan");
        end

        PropDist = [PropDist;CondPropDistMeanSubjs];
        PropAng  = [PropAng;CondPropAngleMeanSubjs];
    
    end

    % remove nans if the row after merging data
    % dataout = removeNanRows(dataout);
    % remove participants for which we were unable to fit the model
    excludedParticipants = removeExcludedParticipants(Glampi_data);

    % In this data format each participant is a column and the
    % environmental conditions are per rows.

    PropDist = PropDist(:,excludedParticipants);
    PropAng  = PropAng(:, excludedParticipants);
end

%% --------------------------------------------------------------------- 
function merged = MergeMCI(MCIUnkData, MCINegData, MCIPosData)
    merged = [MCIUnkData MCINegData MCIPosData];
end

%% ---------------------------------------------------------------------
% Remove participants based on the ability of the GLAMPI model to fit the parameters
function logicalResults = removeExcludedParticipants(groupDataParams)

numConds = 3; % environmental conditions

groupDataParamsAllConds = [];

for trial_filter=1:numConds
    %% extract data
    tempParams = groupDataParams{trial_filter}(:,1); %Can use any of the fitted params
    groupDataParamsAllConds = [groupDataParamsAllConds,tempParams];
end

isnanMatrix = isnan(groupDataParamsAllConds);  % create a logical matrix indicating NaNs
logicalResults = ~all(isnanMatrix, 2);  % find rows with all NaNs

end