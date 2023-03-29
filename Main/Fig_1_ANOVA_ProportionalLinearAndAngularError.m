%% Script to create output for Fig. S8 - parameter comparisons between Young, healthy Elderly, MCI negative
% Zilong Ji, UCL, 2022 zilong.ji@ucl.ac.uk
% Fits the model on all of the groups and run a two-way Anova
% (group*condition) on Young Controls, Elderly Controls and MCI negative
% Output: for each parameter fitted by the model output one boxplot with
% the three groups and performance splitted by environmental condition and
% a boxplot with three groups performance averaged across environmental
% conditions

% Preparing the data
VAM_PrepareBaseConfig;

% Preprocessing the data
VAM_PreprocessData;

% Model fitting
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   100; % skipping the model fitting as we are only interested in the behavioural analysis

VAM;

% Generating color scheme for our paper
ColorPattern;

%% Preparing output
config.ResultFolder = pwd + "/Output/Fig1";
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Collecting information from output
[YoungControlsPropDist, YoungControlsPropAng]     = getProportionalLinearAndAngularError(YoungControls);
[HealthyControlsPropDist, HealthyControlsPropAng] = getProportionalLinearAndAngularError(HealthyControls);
[MCIUnkPropDist, MCIUnkPropAng]                   = getProportionalLinearAndAngularError(MCIUnk);
[MCINegPropDist, MCINegPropAng]                   = getProportionalLinearAndAngularError(MCINeg);
[MCIPosPropDist, MCIPosPropAng]                   = getProportionalLinearAndAngularError(MCIPos);

% TwowayAnova Analysis
config.type = "ProportionalDistance";
[anova_tab_dist,multicomp_tab1_dist,multicomp_tab2_dist, multicomp_tab12_dist]     = TwowayAnovaAllGroups(YoungControlsPropDist, HealthyControlsPropDist, MCIPosPropDist, MCINegPropDist, MCIUnkPropDist, config);

% TwowayAnova Analysis
config.type = "ProportionalAngle";
[anova_tab_angle,multicomp_tab1_angle,multicomp_tab2_angle, multicomp_tab12_angle] = TwowayAnovaAllGroups(YoungControlsPropAng, HealthyControlsPropAng, MCIPosPropAng, MCINegPropAng, MCIUnkPropAng, config);

% Display means and std for each group
disp("%%%%%%%%% Proportional Distance Error %%%%%%%%%");
disp(["Young : " num2str(mean(mean(YoungControlsPropDist,1,"omitnan"),"omitnan")) " +- " std(mean(YoungControlsPropDist,1,"omitnan"),"omitnan")]);
disp(["Elderly : " num2str(mean(mean(HealthyControlsPropDist,1,"omitnan"),"omitnan")) " +- " std(mean(HealthyControlsPropDist,1,"omitnan"),"omitnan")]);
disp(["MCINeg : " num2str(mean(mean(MCIUnkPropDist,1,"omitnan"),"omitnan")) " +- " std(mean(MCIUnkPropDist,1,"omitnan"),"omitnan")]);
disp(["MCIUnk : " num2str(mean(mean(MCINegPropDist,1,"omitnan"),"omitnan")) " +- " std(mean(MCINegPropDist,1,"omitnan"),"omitnan")]);
disp(["MCIPos : " num2str(mean(mean(MCIPosPropDist,1,"omitnan"),"omitnan")) " +- " std(mean(MCIPosPropDist,1,"omitnan"),"omitnan")]);

disp("%%%%%%%%% Proportional Angular Error %%%%%%%%%");
disp(["Young : " num2str(mean(mean(YoungControlsPropAng,1,"omitnan"),"omitnan")) " +- " std(mean(YoungControlsPropAng,1,"omitnan"),"omitnan")]);
disp(["Elderly : " num2str(mean(mean(HealthyControlsPropAng,1,"omitnan"),"omitnan")) " +- " std(mean(HealthyControlsPropAng,1,"omitnan"),"omitnan")]);
disp(["MCINeg : " num2str(mean(mean(MCIUnkPropAng,1,"omitnan"),"omitnan")) " +- " std(mean(MCIUnkPropAng,1,"omitnan"),"omitnan")]);
disp(["MCIUnk : " num2str(mean(mean(MCINegPropAng,1,"omitnan"),"omitnan")) " +- " std(mean(MCINegPropAng,1,"omitnan"),"omitnan")]);
disp(["MCIPos : " num2str(mean(mean(MCIPosPropAng,1,"omitnan"),"omitnan")) " +- " std(mean(MCIPosPropAng,1,"omitnan"),"omitnan")]);

% Proportional Distance Error
plotInfo.defaultTextSize = 20;
plotInfo.defaultLineSize = 2;
plotInfo.barFaceAlpha = 0.5;
plotInfo.barLineWidth = 1.5;
plotInfo.scatterFaceAlpha = 0.2;
plotInfo.scatterEdgeAlpha = 0.8;
plotInfo.scatterDataSize = 32;
plotInfo.errorBarWidth = 2.0;
plotInfo.type = "ProportionalDistance";
plotInfo.YLabel = "Actual distance / correct distance";
plotInfo.visible = "off";

% These values have been set after manually checking the output of the two way anova
% This adds only stars on top of the groups
% We are not going to set sigstar here as we are going to group them to
% make visualization more easy
plotInfo.addStar = [nan nan nan nan nan];
plotInfo.starXValues = [1 2 3 4 5];
plotInfo.starYValues = [nan nan nan nan 1.5];
plotInfo.starTextSize = [35 35 35 35 35];
% This adds sigstar bars using adjustablesigstar (if set to 1)
plotInfo.addSigmaStarbars = 0;
plotInfo.sigmaStarBars = [4 5; 3 5; 2 5; 1 5; 1 2];
plotInfo.sigmaStarBarsPValues = [0.001; 0.001; 0.001; 0.001; 0.02];
plotInfo.sigmaStarLineWidth = 2.5;
plotInfo.sigmaStarTextSize  = 20;
plotInfo.sigmaBarSeparation = 0.075;
plotInfo.yLim = [0 1.75];
plotInfo.ticksStep = 0.5;

plotInfo.yLim = [0 1.75];
plotInfo.ticksStep = 0.5;

plotMergedBarScatter(mean(YoungControlsPropDist,1,'omitnan'), mean(HealthyControlsPropDist,1,'omitnan'), mean(MCIPosPropDist,1,'omitnan'), mean(MCINegPropDist,1,'omitnan'), mean(MCIUnkPropDist,1,'omitnan'), anova_tab_dist, multicomp_tab1_dist, config, plotInfo);

% Proportional Angular Error
plotInfo.defaultTextSize = 20;
plotInfo.defaultLineSize = 2;
plotInfo.barFaceAlpha = 0.5;
plotInfo.barLineWidth = 1.5;
plotInfo.scatterFaceAlpha = 0.2;
plotInfo.scatterEdgeAlpha = 0.8;
plotInfo.scatterDataSize = 32;
plotInfo.errorBarWidth = 2.0;
plotInfo.type = "ProportionalAngle";
plotInfo.YLabel = "Actual angle / correct angle";
plotInfo.visible = "off";
% These values have been set after checking the output of the two way anova
% This adds only stars on top of the groups
plotInfo.addStar = [nan nan nan nan nan];
plotInfo.starXValues = [1 2 3 4 5];
plotInfo.starYValues = [nan nan nan nan 2.1];
plotInfo.starTextSize = [35 35 35 35 35];
% This adds sigstar bars using adjustablesigstar (if set to 1)
plotInfo.addSigmaStarbars = 0;
plotInfo.sigmaStarBars = [4 5; 3 5; 2 5; 1 5];
plotInfo.sigmaStarBarsPValues = [0.001; 0.001; 0.001; 0.001];
plotInfo.sigmaStarLineWidth = 2.5;
plotInfo.sigmaStarTextSize  = 20;
plotInfo.sigmaBarSeparation = 0.06;
plotInfo.yLim = [0 2.5];
plotInfo.ticksStep = 0.5;

plotInfo.yLim = [0 2.5];
plotMergedBarScatter(mean(YoungControlsPropAng,1,'omitnan'), mean(HealthyControlsPropAng,1,'omitnan'), mean(MCIPosPropAng,1,'omitnan'), mean(MCINegPropAng,1,'omitnan'), mean(MCIUnkPropAng,1,'omitnan'), anova_tab_angle, multicomp_tab1_angle, config, plotInfo);

% Final cleanup to leave workspace as the end of the Preprocessing stage.
% Remove if you want to take a look at the output data.
clearvars -except config YoungControls HealthyControls MCINeg MCIPos MCIUnk anova_tab_angle anova_tab_dist multicomp_tab1_dist multicomp_tab2_dist multicomp_tab12_dist multicomp_tab1_angle multicomp_tab2_angle multicomp_tab12_angle

%% --------------------------------------------------------------------- 
function [PropDist, PropAng] = getProportionalLinearAndAngularError(Group)
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
end

%% ---------------------------------------------------------------------
function plotMergedBarScatter(YoungData, HealthyOldData, MCIPosData, MCINegData, MCIUnkData, ANOVA_tab, ANOVA_multcomp_group, config, plotInfo)
% plot Bar scatter of data provided
% Data has to be formatted to be one line for each single group. Group
% plotting is always young, older, mci unk, mci neg, mci pos

% Recombining the data
AllData{1} = YoungData; AllData{2} = HealthyOldData; AllData{3} = MCIUnkData; AllData{4} = MCINegData; AllData{5} = MCIPosData;

mean_all   = [mean(YoungData, 2,"omitnan")';%avearge over the row;
    mean(HealthyOldData, 2,"omitnan")';
    mean(MCIUnkData, 2,"omitnan")';
    mean(MCINegData, 2,"omitnan")';
    mean(MCIPosData, 2,"omitnan")'];

std_all = [std(YoungData,0,2,"omitnan")'      ./sqrt(sum(~isnan(YoungData),2,"omitnan")');
    std(HealthyOldData,0,2,"omitnan")' ./sqrt(sum(~isnan(HealthyOldData),2,"omitnan")');
    std(MCIUnkData,0,2,"omitnan")'     ./sqrt(sum(~isnan(MCIUnkData),2,"omitnan")');
    std(MCINegData,0,2,"omitnan")'     ./sqrt(sum(~isnan(MCINegData),2,"omitnan")');
    std(MCIPosData,0,2,"omitnan")'     ./sqrt(sum(~isnan(MCIPosData),2,"omitnan")')];

f = figure('visible',plotInfo.visible,'Position', [100 100 850 500]);

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',plotInfo.defaultTextSize)
set(0,'DefaultTextFontSize',plotInfo.defaultTextSize)

%% Bar Plotting
b = bar(mean_all, 'grouped');
b.FaceColor = 'flat';
b.FaceAlpha = plotInfo.barFaceAlpha;
b.LineWidth = plotInfo.barLineWidth;
for j=1:length(mean_all)
    b.CData(j,:) = config.color_scheme_group(j,:);
end

hold on

x = b.XEndPoints;

scatter_facealpha=plotInfo.scatterFaceAlpha;
scatter_edgealpha=plotInfo.scatterEdgeAlpha;
scatter_datasize=plotInfo.scatterDataSize;

for j=1:length(mean_all)
    % Scatter data
    scatter(x(j), AllData{j},30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor',config.color_scheme_group(j,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha, SizeData=scatter_datasize);
end

hold on

% Error Bars
errorbar(x',mean_all,std_all,'k','LineStyle','None', 'LineWidth', plotInfo.errorBarWidth, 'CapSize', 14);

hold on
% Reference Line
yline(1, 'LineStyle','-.', 'LineWidth', 1.5, 'Color', config.color_scheme_npg(1,:), 'Alpha', 0.8);

% Figure post-processing
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...
    'XTick'       , (1:5),...
    'XTickLabel'  , ["Young", "Elderly", "MCI Unk", "MCI Neg", "MCI Pos"],...
    'LineWidth'   , .5        );

sigstaroptions.textSize      = plotInfo.sigmaStarTextSize;
sigstaroptions.lineWidth     = plotInfo.sigmaStarLineWidth;
sigstaroptions.barSeparation = plotInfo.sigmaBarSeparation;

if(plotInfo.addSigmaStarbars == 1)
    for j = 1:height(plotInfo.sigmaStarBars)
        adjustablesigstar(plotInfo.sigmaStarBars(j,:),plotInfo.sigmaStarBarsPValues(j),0,sigstaroptions);
    end
end

ylabel(plotInfo.YLabel);
ylim(plotInfo.yLim);
yticks(0:plotInfo.ticksStep:plotInfo.yLim(2));

xtickangle(35)

ax = gca;
ax.LineWidth = plotInfo.defaultLineSize;

hold off;

% export figure
exportgraphics(f,config.ResultFolder+"/BarMergedNormedReturn" + plotInfo.type + ".png",'Resolution',300);
exportgraphics(f,config.ResultFolder+"/BarMergedNormedReturn" + plotInfo.type + ".pdf",'Resolution',300, 'ContentType','vector');

end