%% Script to create output for Supp Fig. 2 - analysis of time for completing outbound duration
% Andrea Castegnaro, UCL, 2022 andrea.castegnaro@ucl.ac.uk
% Display outbound path duration for different groups.
% Calculate difference between group means. For details of how duration has
% been calculated please refer to CalculateTrackingPath function.

% Preparing the data
GLAMPI_PrepareBaseConfig

% Preprocessing the data
GLAMPI_PreprocessData

% Model fitting
rng("default");
% No need to filter between different environmental condition
config.useTrialFilter = false;
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   100; % skip model fitting

GLAMPI

% Preparing output
config.ResultFolder = pwd + "/Output/S2/TimeAnalysis";
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Generating color scheme for our paper
ColorPattern;

%% Collecting information from output
YoungOutboundDuration = extractTimeData(YoungControls);
HealthyControlsOutboundDuration = extractTimeData(HealthyControls);
MCINegOutboundDuration = extractTimeData(MCINeg);
MCIPosOutboundDuration = extractTimeData(MCIPos);
MCIUnkOutboundDuration = extractTimeData(MCIUnk);

% Analysis and plotting
close all;

% Parameters set for controlling visual output
plotInfo.defaultLineSize = 1.7;
plotInfo.titleFontSize = 12;
plotInfo.labelSize = 12;
plotInfo.axisSize = 10;
plotInfo.MarkerSize = 20;
plotInfo.MarkerAlpha = 0.5;
plotInfo.PatchAlpha = 0.7;
plotInfo.yLim = [8 20];
plotInfo.xLim = [0.5 5.5];
plotInfo.medianColor = [0.4 0.4 0.4];
plotInfo.medianWidth = 1.3;
plotInfo.meanMarkerSize = 30;
plotInfo.sigmaStarLineWidth = 2.5;
plotInfo.sigmaStarTextSize  = 20;
plotInfo.sigmaBarSeparation = 0.04;
plotInfo.FigurePosition = [200 200 280 250];

% Generate example data
max_element = max([numel(YoungOutboundDuration),...
    numel(HealthyControlsOutboundDuration),...
    numel(MCINegOutboundDuration),...
    numel(MCIPosOutboundDuration),...
    numel(MCIUnkOutboundDuration)]);
YoungOutboundDuration(end+1 : max_element) = nan;
HealthyControlsOutboundDuration(end+1 : max_element) = nan;
MCIUnkOutboundDuration(end+1 : max_element) = nan;
MCINegOutboundDuration(end+1 : max_element) = nan;
MCIPosOutboundDuration(end+1 : max_element) = nan;

OutboundData = [YoungOutboundDuration HealthyControlsOutboundDuration MCIUnkOutboundDuration MCINegOutboundDuration MCIPosOutboundDuration];

DataMeans = mean(OutboundData,1,"omitnan");
DataSems = std(OutboundData,1,"omitnan")./sqrt(sum(~isnan(OutboundData),1));

xDatameans = 1:length(DataMeans);

% Calculating anova
OutboundDataAnova = [YoungOutboundDuration; HealthyControlsOutboundDuration; MCIUnkOutboundDuration; MCINegOutboundDuration; MCIPosOutboundDuration];

OutboundDataAnovaGroups = [repmat({'Young'}, size(YoungOutboundDuration, 1), 1); ...
          repmat({'Healthy Older'}, size(HealthyControlsOutboundDuration, 1), 1); ...
          repmat({'MCI Unknown'}, size(MCIUnkOutboundDuration, 1), 1); ...
          repmat({'MCI Negative'}, size(MCINegOutboundDuration, 1), 1); ...
          repmat({'MCI Positive'}, size(MCIPosOutboundDuration, 1), 1)];

% Anova
[p, tbl, stats] = anova1(OutboundDataAnova, OutboundDataAnovaGroups, 'off');
disp("Anova on Time - 1 Young 2 Elderly 3 Mci unk 4 Mci neg 5 Mci pos");
tbl
% Multiple comparisons
[results, means, ~, ~] = multcompare(stats, "Alpha", 0.01, "CType","bonferroni", 'Display','off');

% Reporting multiple comparisons results
disp("Anova multiple comparisons results - 1 Young 2 Elderly 3 Mci unk 4 Mci neg 5 Mci pos");
results

currFig = figure("Position",plotInfo.FigurePosition,"Visible","off");

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

hold on;

% Create box plot
boxplots = boxplot(OutboundData,"Notch","on","symbol","","Colors",config.color_scheme_group, "Widths",0.6);
boxes = findobj(gca,'Tag','Box');
boxes = flip(boxes);
medians = findobj(gca,'Tag','Median');
medians = flip(medians);

for j = 1:length(boxes)
    pp = patch(get(boxes(j),'XData'), get(boxes(j), 'YData'), config.color_scheme_group(j,:), 'FaceAlpha', plotInfo.PatchAlpha);
    pp.LineStyle = "none";
end
clear j;

for j = 1:length(medians)
    plot(medians(j).XData,medians(j).YData,Color=plotInfo.medianColor, LineStyle="-", LineWidth=plotInfo.medianWidth);
end
clear j;

for j = 1:width(OutboundData)
    sh = scatter(j*ones(height(OutboundData),1), OutboundData(:,j));
    sh.SizeData = plotInfo.MarkerSize;
    sh.MarkerEdgeColor = "none";
    sh.MarkerFaceColor = config.color_scheme_group(j,:);
    sh.MarkerFaceAlpha = plotInfo.MarkerAlpha;
end
clear j;

ax_errorBar = errorbar(xDatameans,DataMeans,DataSems);
ax_errorBar.Color = [0 0 0];
ax_errorBar.LineWidth = 3;
ax_errorBar.LineStyle = "none";

sc_means = scatter(xDatameans,DataMeans);
sc_means.Marker = "diamond";
sc_means.SizeData = plotInfo.meanMarkerSize;
sc_means.MarkerFaceAlpha = 1;
sc_means.MarkerFaceColor = "white";
sc_means.MarkerEdgeColor = "none";

sigstaroptions.textSize      = plotInfo.sigmaStarTextSize;
sigstaroptions.lineWidth     = plotInfo.sigmaStarLineWidth;
sigstaroptions.barSeparation = plotInfo.sigmaBarSeparation;

adjustablesigstar([1 2],0.001,0,sigstaroptions);
adjustablesigstar([2 3],0.05,0,sigstaroptions);

hold off;

ylim(plotInfo.yLim);
xlim(plotInfo.xLim);

ylabel("Duration (s)");

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1]);


ax = gca;
ax.XAxis.TickLabels = {'Young', 'Elderly', 'MCI Unk', 'MCI Neg', 'MCI Pos'};

ax.LineWidth = plotInfo.defaultLineSize;
ax.XLabel.FontSize = plotInfo.labelSize;
ax.YLabel.FontSize = plotInfo.labelSize;
ax.XAxis.FontSize = plotInfo.axisSize;
ax.YAxis.FontSize = plotInfo.axisSize;

exportgraphics(currFig,config.ResultFolder+"/OutboundPathDuration.png",'Resolution',300);
exportgraphics(currFig,config.ResultFolder+"/OutboundPathDuration.pdf",'Resolution',300, 'ContentType','vector');

% Final cleanup to leave workspace as the end of the Preprocessing stage.
% Remove if you want to take a look at the output data.
clearvars -except config YoungControls HealthyControls MCINeg MCIPos MCIUnk

%% ------------------------------------------------------------------
function dataout = extractTimeData(groupData)

    dataout = [];
    
    for i_participant = 1:length(groupData.Results.L1Dur)
        if(~isempty(groupData.Results.L1Dur{i_participant}))
            l1_duration = cell2mat(groupData.Results.L1Dur{i_participant});
            l2_duration = cell2mat(groupData.Results.L2Dur{i_participant});
            dataout = [dataout; mean(l1_duration+l2_duration,'omitnan')];
        end
    end
end