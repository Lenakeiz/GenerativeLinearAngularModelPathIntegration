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
config.ResultFolder = pwd + "/Output/Fig1/YoungvsHealthyOldvsMCIMerged_Cut9trial";
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Collecting information from output
[YoungControlsPropDist, YoungControlsPropAng]     = getProportionalLinearAndAngularError(YoungControls, YoungControls.Results.estimatedParams, true);
[HealthyControlsPropDist, HealthyControlsPropAng] = getProportionalLinearAndAngularError(HealthyControls, HealthyControls.Results.estimatedParams, false);
[MCIUnkPropDist, MCIUnkPropAng]                   = getProportionalLinearAndAngularError(MCIUnk, MCIUnk.Results.estimatedParams, false);
[MCINegPropDist, MCINegPropAng]                   = getProportionalLinearAndAngularError(MCINeg, MCINeg.Results.estimatedParams, false);
[MCIPosPropDist, MCIPosPropAng]                   = getProportionalLinearAndAngularError(MCIPos, MCIPos.Results.estimatedParams, false);

% MergeMCI
MCIMergedPropDist = MergeMCI(MCIUnkPropDist, MCINegPropDist, MCIPosPropDist);
MCIMergedPropAng  = MergeMCI(MCIUnkPropAng, MCINegPropAng, MCIPosPropAng);

% TwowayAnova Analysis
config.type = "ProportionalDistance";
[anova_tab_dist,multicomp_tab1_dist,multicomp_tab2_dist, multicomp_tab12_dist]     = TwowayAnova_Behavioural_YoungvsElderlyvsMCI(YoungControlsPropDist, HealthyControlsPropDist, MCIMergedPropDist, config);

% TwowayAnova Analysis
config.type = "ProportionalAngle";
[anova_tab_angle,multicomp_tab1_angle,multicomp_tab2_angle, multicomp_tab12_angle] = TwowayAnova_Behavioural_YoungvsElderlyvsMCI(YoungControlsPropAng, HealthyControlsPropAng, MCIMergedPropAng, config);

% Display means and std for each group
disp("%%%%%%%%% Proportional Distance Error %%%%%%%%%");
disp(["Young : " num2str(mean(mean(YoungControlsPropDist,1,"omitnan"),"omitnan")) " +- " std(mean(YoungControlsPropDist,1,"omitnan"),"omitnan")]);
disp(["Elderly : " num2str(mean(mean(HealthyControlsPropDist,1,"omitnan"),"omitnan")) " +- " std(mean(HealthyControlsPropDist,1,"omitnan"),"omitnan")]);
disp(["MCIMerged : " num2str(mean(mean(MCIMergedPropDist,1,"omitnan"),"omitnan")) " +- " std(mean(MCIMergedPropDist,1,"omitnan"),"omitnan")]);

TTest_NominalValue(mean(YoungControlsPropDist,1,"omitnan")',1,'Young', 'both');
TTest_NominalValue(mean(HealthyControlsPropDist,1,"omitnan")',1,'Elderly', 'both');
TTest_NominalValue(mean(MCIMergedPropDist,1,"omitnan")',1,'MCI', 'both');

disp("%%%%%%%%% Proportional Angular Error %%%%%%%%%");
disp(["Young : " num2str(mean(mean(YoungControlsPropAng,1,"omitnan"),"omitnan")) " +- " std(mean(YoungControlsPropAng,1,"omitnan"),"omitnan")]);
disp(["Elderly : " num2str(mean(mean(HealthyControlsPropAng,1,"omitnan"),"omitnan")) " +- " std(mean(HealthyControlsPropAng,1,"omitnan"),"omitnan")]);
disp(["MCIMerged : " num2str(mean(mean(MCIMergedPropAng,1,"omitnan"),"omitnan")) " +- " std(mean(MCIMergedPropAng,1,"omitnan"),"omitnan")]);

TTest_NominalValue(mean(YoungControlsPropAng,1,"omitnan")',1,'Young', 'both');
TTest_NominalValue(mean(HealthyControlsPropAng,1,"omitnan")',1,'Elderly', 'both');
TTest_NominalValue(mean(MCIMergedPropAng,1,"omitnan")',1,'MCI', 'both');

% Proportional Distance Error
plotInfo.defaultTextSize = 20;
plotInfo.defaultLineSize = 2;
plotInfo.barFaceAlpha = 0.5;
plotInfo.barLineWidth = 0.3;
plotInfo.barWidth = 0.5;
plotInfo.scatterFaceAlpha = 0.5;
plotInfo.scatterEdgeAlpha = 0.8;
plotInfo.scatterDataSize = 32;
plotInfo.errorBarWidth = 2.0;
plotInfo.visible = "off";
plotInfo.dimensions = [100 100 250 250];

plotInfo.type = "ProportionalDistance";
plotInfo.YLabel = "Actual distance / correct distance";
plotInfo.yLim = [0 1.75];
plotInfo.ticksStep = 0.5;

BoxPlotOfFittedParamMergeCondition(YoungControlsPropDist, HealthyControlsPropDist, MCIMergedPropDist, anova_tab_dist, multicomp_tab1_dist, config, plotInfo);

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
BoxPlotOfFittedParamMergeCondition(YoungControlsPropAng, HealthyControlsPropAng, MCIMergedPropAng, anova_tab_angle, multicomp_tab1_angle, config, plotInfo);

% Final cleanup to leave workspace as the end of the Preprocessing stage.
% Remove if you want to take a look at the output data.
% clearvars -except config YoungControls HealthyControls MCINeg MCIPos MCIUnk anova_tab_angle anova_tab_dist multicomp_tab1_dist multicomp_tab2_dist multicomp_tab12_dist multicomp_tab1_angle multicomp_tab2_angle multicomp_tab12_angle

%% ---------------------------------------------------------------------
function BoxPlotOfFittedParamMergeCondition(YoungData, HealthyData, MCIData, ANOVA_tab, ANOVA_multcomp_group, config, plotInfo)

YoungParamMean      = mean(YoungData, 1, "omitnan")';
HealthyOldParamMean = mean(HealthyData, 1, "omitnan")';
MCIParamMean        = mean(MCIData, 1, "omitnan")';

f = figure('visible', plotInfo.visible,'Position', plotInfo.dimensions);

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12)

colorForYoung = config.color_scheme_npg(3,:);
colorForElderly = config.color_scheme_npg(5,:);
colorForMCI = config.color_scheme_npg(2,:);

% parameters set for controlling visual output
whisker_value               =   1.5;
box_lineWidth               =   plotInfo.barLineWidth;
box_widths_value            =   plotInfo.barWidth;
box_color_transparency      =   plotInfo.barFaceAlpha; %faceAlpha
median_lineWidth            =   2;
median_color                =   'k';
scatter_jitter_value        =   0.3;
scatter_markerSize          =   plotInfo.scatterDataSize;
scatter_marker_edgeColor    =   'k';
scatter_marker_edgeWidth    =   0.5;
scatter_color_transparency  =   plotInfo.scatterFaceAlpha; %faceAlpha
mean_scatter_multiplier     =   2.0;

hold on
%% Boxplot for each column in Young
bp1 = boxplot(YoungParamMean, ...
    'Whisker',whisker_value, ...
    'symbol','', ... %symbol ='' making outlier invisible
    'Color','k', ...
    'Notch','on', ...
    'widths',box_widths_value,...
    'positions', 1);
set(bp1,'linewidth',box_lineWidth);

hold on
%% Boxplot for each column in Healthy Elderly
bp2 = boxplot(HealthyOldParamMean, ...
    'Whisker',whisker_value, ...
    'symbol','', ... %symbol ='' making outlier invisible
    'Color','k', ...
    'Notch','on', ...
    'widths',box_widths_value,...
    'positions', 2);
set(bp2,'linewidth',box_lineWidth);

hold on
%% Boxplot for each column in MCI negative
bp3 = boxplot(MCIParamMean, ...
    'Whisker',whisker_value, ...
    'symbol','', ... %symbol ='' making outlier invisible
    'Color','k', ...
    'Notch','on', ...
    'widths',box_widths_value,...
    'positions', 3);
set(bp3,'linewidth',box_lineWidth);

%% Boxplot visual changes
% matlab has a lifo system for these
h = findobj(gca,'Tag','Box');
%get the Young box
patch(get(h(3),'XData'),get(h(3),'YData'),colorForYoung,'FaceAlpha',box_color_transparency);
%get the HelthyOld box
patch(get(h(2),'XData'),get(h(2),'YData'),colorForElderly,'FaceAlpha',box_color_transparency);
%get the MCI box
patch(get(h(1),'XData'),get(h(1),'YData'),colorForMCI,'FaceAlpha',box_color_transparency);

%% Median visual change
h=findobj(gca,'tag','Median');
for i = 1:length(h)
    h(i).LineWidth = median_lineWidth;
    h(i).Color = median_color;
end

%% Scatter plot for data and mean (Young)
num_points = size(YoungParamMean,1);
hold on
x = 1*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
scatter(x, YoungParamMean, scatter_markerSize, ...
    'filled', ...
    'o', ...
    'MarkerEdgeColor',scatter_marker_edgeColor, ...
    'MarkerFaceColor',colorForYoung, ...
    'MarkerFaceAlpha',scatter_color_transparency,...
    'LineWidth',scatter_marker_edgeWidth);

hold on
%add errorbar
mean_Young = mean(YoungParamMean, "omitnan");
sem_Young = std(YoungParamMean,"omitnan")./sqrt(length(YoungParamMean(~isnan(YoungParamMean))));
errorbar(1,mean_Young,sem_Young,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18);
hold on
%add mean point
scatter(1, mean_Young, mean_scatter_multiplier*scatter_markerSize, 'd',...
    'filled','MarkerEdgeColor','k', ...
    'MarkerFaceColor','w', ...
    'LineWidth',scatter_marker_edgeWidth);

%% Scatter plot for data and mean (Healthy Elderly)
num_points = length(HealthyOldParamMean);
hold on
x = 2*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
scatter(x, HealthyOldParamMean, scatter_markerSize, ...
    'filled', ...
    'o', ...
    'MarkerEdgeColor',scatter_marker_edgeColor, ...
    'MarkerFaceColor',colorForElderly, ...
    'MarkerFaceAlpha',scatter_color_transparency,...
    'LineWidth',scatter_marker_edgeWidth);

%add errorbar
mean_Hold = mean(HealthyOldParamMean, "omitnan");
sem_Hold = std(HealthyOldParamMean, "omitnan")./sqrt(length(HealthyOldParamMean(~isnan(HealthyOldParamMean))));
errorbar(2,mean_Hold,sem_Hold,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18);
hold on
%add mean point
scatter(2, mean_Hold, mean_scatter_multiplier*scatter_markerSize, 'd',...
    'filled','MarkerEdgeColor','k', ...
    'MarkerFaceColor','w', ...
    'LineWidth',scatter_marker_edgeWidth);

%% Scatter plot for data and mean (MCI unknown)
num_points = length(MCIParamMean);
hold on
x = 3*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
scatter(x, MCIParamMean, scatter_markerSize, ...
    'filled', ...
    'o', ...
    'MarkerEdgeColor',scatter_marker_edgeColor, ...
    'MarkerFaceColor',colorForMCI, ...
    'MarkerFaceAlpha',scatter_color_transparency,...
    'LineWidth',scatter_marker_edgeWidth);

%add mean + errorbar
mean_MCI = mean(MCIParamMean, "omitnan");
sem_MCI = std(MCIParamMean, "omitnan")./sqrt(length(MCIParamMean(~isnan(MCIParamMean))));
errorbar(3,mean_MCI,sem_MCI,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18);
hold on
%add mean point
scatter(3, mean_MCI, mean_scatter_multiplier*scatter_markerSize, 'd',...
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
    'XTick'       , (1:3),...
    'XLim'        , [0.5, 3.5],...
    'XTickLabel'  , {'Young','Elderly','MCI'},...
    'LineWidth'   , 1.0        );

ylabel(plotInfo.YLabel);
ylim(plotInfo.yLim);
yticks(0:plotInfo.ticksStep:plotInfo.yLim(2));

% Extract pvalues from multiple comparison of group effect. Adding
% it to the figure
multicomp_result = ANOVA_multcomp_group;
% 1       2             3      
% Young   HealthyOld    MCI
PvalueYoungvsHealthyOld = multicomp_result(1,6); % Young vs. HealthyOld
PvalueHealthyOldvsMCI = multicomp_result(3,6); % HealthyOld v.s. MCI unk
PvalueYoungvsMCI = multicomp_result(2,6);


%% Add significance bars
AllP = [PvalueYoungvsMCI,PvalueYoungvsHealthyOld,PvalueHealthyOldvsMCI];
Xval = [[1,3];[1,2];[2,3]];
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
function [PropDist, PropAng] = getProportionalLinearAndAngularError(Group, Glampi_data, cuttrial)
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
            if cuttrial & ~isempty(CondPropDist{1,subj})
                cpd = cell2mat(CondPropDist{1,subj});
                cpd = cpd(1:9);
                cpa = cell2mat(CondPropAng{1,subj});
                cpa = cpa(1:9);
                CondPropDistMeanSubjs(1,subj)  = mean(cpd,"omitnan");
                CondPropAngleMeanSubjs(1,subj) = mean(cpa,"omitnan");
            else
                CondPropDistMeanSubjs(1,subj)  = mean(cell2mat(CondPropDist{1,subj}),"omitnan");
                CondPropAngleMeanSubjs(1,subj) = mean(cell2mat(CondPropAng{1,subj}),"omitnan");
            end
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