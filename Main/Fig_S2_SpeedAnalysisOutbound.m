%% Script to create output for Supp Fig. 2 - analysis of walking speed during outbound path
% Andrea Castegnaro, UCL, 2022 andrea.castegnaro@ucl.ac.uk
% Display walking speed for different groups during the outbound path.
% Calculate difference between group means. For details of speed
% calculation please refer to CalculateTrackingPath function.

% Preparing the data
VAM_PrepareBaseConfig

% Preprocessing the data
VAM_PreprocessData

% Model fitting
rng("default");
% No need to filter between different environmental condition
config.useTrialFilter = false;
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   100; % skip model fitting

VAM

% Preparing output
config.ResultFolder = pwd + "/Output/S2/SpeedAnalysis";
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Generating color scheme for our paper
ColorPattern;

%% Analysis and plotting
% Parameters set for controlling visual output
plotInfo.defaultLineSize = 1.7;
plotInfo.titleFontSize = 12;
plotInfo.labelSize = 12;
plotInfo.axisSize = 10;
plotInfo.MarkerSize = 20;
plotInfo.MarkerAlpha = 0.5;
plotInfo.PatchAlpha = 0.7;
plotInfo.yLim = [0.2 1.0];
plotInfo.xLim = [0.5 5.5];
plotInfo.medianColor = [0.4 0.4 0.4];
plotInfo.medianWidth = 1.3;
plotInfo.meanMarkerSize = 30;
plotInfo.sigmaStarLineWidth = 2.5;
plotInfo.sigmaStarTextSize  = 20;
plotInfo.sigmaBarSeparation = 0.06;
plotInfo.FigurePosition = [200 200 280 250];

% Averaging across segments L1 and L2 of the outbound path
YoungSpeed = extractSpeedData(YoungControls);
HealthyControlsSpeed = extractSpeedData(HealthyControls);
MCINegSpeed = extractSpeedData(MCINeg);
MCIPosSpeed = extractSpeedData(MCIPos);
MCIUnkSpeed = extractSpeedData(MCIUnk);

% Generate example data
max_element = max([numel(YoungSpeed), numel(HealthyControlsSpeed), numel(MCINegSpeed), numel(MCIPosSpeed), numel(MCIUnkSpeed)]);
YoungSpeed(end+1 : max_element) = nan;
HealthyControlsSpeed(end+1 : max_element) = nan;
MCIUnkSpeed(end+1 : max_element) = nan;
MCINegSpeed(end+1 : max_element) = nan;
MCIPosSpeed(end+1 : max_element) = nan;

OutboundData = [YoungSpeed HealthyControlsSpeed MCIUnkSpeed MCINegSpeed MCIPosSpeed];

DataMeans = mean(OutboundData,1,"omitnan");
DataSems = std(OutboundData,1,"omitnan")./sqrt(sum(~isnan(OutboundData),1));

xDatameans = 1:length(DataMeans);

%
OutboundDataAnova = [YoungSpeed; HealthyControlsSpeed; MCIUnkSpeed; MCINegSpeed; MCIPosSpeed];

OutboundDataAnovaGroups = [repmat({'Young'}, size(YoungSpeed, 1), 1); ...
          repmat({'Healthy Older'}, size(HealthyControlsSpeed, 1), 1); ...
          repmat({'MCI Unknown'}, size(MCIUnkSpeed, 1), 1); ...
          repmat({'MCI Negative'}, size(MCINegSpeed, 1), 1); ...
          repmat({'MCI Positive'}, size(MCIPosSpeed, 1), 1)];

[p, tbl, stats] = anova1(OutboundDataAnova, OutboundDataAnovaGroups, 'off');
[results, means, ~, ~] = multcompare(stats, 'alpha', 0.05,'Display','off');

% close all
close all;

currFig = figure("Position",plotInfo.FigurePosition, Visible="off");

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

adjustablesigstar([1 2],0.05,0,sigstaroptions);

hold off;

ylim(plotInfo.yLim);
xlim(plotInfo.xLim);

ylabel("Speed (m/s)");

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

exportgraphics(currFig,config.ResultFolder+"/OutboundPathSpeed.png",'Resolution',300);
exportgraphics(currFig,config.ResultFolder+"/OutboundPathSpeed.pdf",'Resolution',300, 'ContentType','vector');

% Final cleanup to leave workspace as the end of the Preprocessing stage.
% Remove if you want to take a look at the output data.
clearvars -except config YoungControls HealthyControls MCINeg MCIPos MCIUnk

%% ------------------------------------------------------------------
function dataout = extractSpeedData(groupData)

    dataout = nan(length(groupData.TrackedL1),2);
    
    for i_participant = 1:length(groupData.TrackedL1)
        % Extracting data calculated from the first segment
        speedTemp = nan(height(groupData.TrackedL1{i_participant}),1);
        for i_trial = 1:height(groupData.TrackedL1{i_participant})
            speedTemp(i_trial) = mean(groupData.TrackedL1{i_participant}{i_trial}.Smoothed_Vel_proj(groupData.TrackedL1{i_participant}{i_trial}.Filtered_Vel_proj),'omitnan');
        end
        dataout(i_participant,1) = mean(speedTemp,'omitnan');

        speedTemp = nan(height(groupData.TrackedL2{i_participant}),1);
        % Extracting data calculated from the second segment
        for i_trial = 1:height(groupData.TrackedL2{i_participant})
            speedTemp(i_trial) = mean(groupData.TrackedL2{i_participant}{i_trial}.Smoothed_Vel_proj(groupData.TrackedL2{i_participant}{i_trial}.Filtered_Vel_proj),'omitnan');
        end
        dataout(i_participant,2) = mean(speedTemp,'omitnan');
    end
    % Calculating average across rows to calculate average speed of participants
    % during the two segments.
    dataout = mean(dataout,2,'omitnan');
end