%% Script to create output for Supp Fig. 2 - analysis of pitch rotation during turn between outbound segments
% Andrea Castegnaro, UCL, 2022 andrea.castegnaro@ucl.ac.uk
% Calculate pitch rotation during the rotation at the outbound segment. A
% different pitch rotation might indicate a difference in angular velocity
% encoding (please see Choi et al., 2023) for details

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
config.ResultFolder = pwd + "/Output/S2/YawAnalysis";
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Generating color scheme for our paper
ColorPattern;

%% Plotting the calculation of the yaw points for some participant
% Unity is y-up world coordinate
drawHeadingPoints(YoungControls.Path{1,1}{1,1}(:,[5 7 6]),YoungControls.Path{1,1}{1,1}(:,1),YoungControls.FlagTrigTimes{1,1}{1,2},YoungControls.FlagTrigTimes{1,1}{1,3});
drawHeadingPoints(HealthyControls.Path{1,1}{1,1}(:,[5 7 6]),HealthyControls.Path{1,1}{1,1}(:,1),HealthyControls.FlagTrigTimes{1,1}{1,2},HealthyControls.FlagTrigTimes{1,1}{1,3});
drawHeadingPoints(MCIUnk.Path{1,1}{1,1}(:,[5 7 6]),MCIUnk.Path{1,1}{1,1}(:,1),MCIUnk.FlagTrigTimes{1,1}{1,2},MCIUnk.FlagTrigTimes{1,1}{1,3});
drawHeadingPoints(MCINeg.Path{1,1}{1,1}(:,[5 7 6]),MCINeg.Path{1,1}{1,1}(:,1),MCINeg.FlagTrigTimes{1,1}{1,2},MCINeg.FlagTrigTimes{1,1}{1,3});
drawHeadingPoints(MCIPos.Path{1,1}{1,1}(:,[5 7 6]),MCIPos.Path{1,1}{1,1}(:,1),MCIPos.FlagTrigTimes{1,1}{1,2},MCIPos.FlagTrigTimes{1,1}{1,3});

%% Extracting the Pitch angle
% Calculating the pitch angle of the head during the outbound path between
% cone 2 and cone 3
YoungData = extractPitchData(YoungControls);
HealthyControlsData = extractPitchData(HealthyControls);
MCINegData = extractPitchData(MCINeg);
MCIPosData = extractPitchData(MCIPos);
MCIUnkData = extractPitchData(MCIUnk);

% Young/Healthy Controls/MCI separately
close all;

% Parameters set for controlling visual output
plotInfo.defaultLineSize = 1.7;
plotInfo.titleFontSize = 12;
plotInfo.labelSize = 12;
plotInfo.axisSize = 10;
plotInfo.MarkerSize = 20;
plotInfo.MarkerAlpha = 0.35;
plotInfo.PatchAlpha = 0.7;
plotInfo.yLim = [-30 20];
plotInfo.xLim = [0.5 5.5];
plotInfo.medianColor = [0.4 0.4 0.4];
plotInfo.medianWidth = 1.3;
plotInfo.meanMarkerSize = 30;
plotInfo.sigmaStarLineWidth = 2.5;
plotInfo.sigmaStarTextSize  = 20;
plotInfo.sigmaBarSeparation = 0.04;
plotInfo.FigurePosition = [200 200 280 250];

% Filling with nans smaller columns
max_element = max([numel(YoungData),...
    numel(HealthyControlsData),...
    numel(MCINegData),...
    numel(MCIPosData),...
    numel(MCIUnkData)]);
YoungData(end+1 : max_element) = nan;
HealthyControlsData(end+1 : max_element) = nan;
MCIUnkData(end+1 : max_element) = nan;
MCINegData(end+1 : max_element) = nan;
MCIPosData(end+1 : max_element) = nan;

OutboundData = [YoungData HealthyControlsData MCIUnkData MCINegData MCIPosData];

DataMeans = mean(OutboundData,1,"omitnan");
DataSems = std(OutboundData,1,"omitnan")./sqrt(sum(~isnan(OutboundData),1));

xDatameans = 1:length(DataMeans);

%
OutboundDataAnova = [YoungData; HealthyControlsData; MCIUnkData; MCINegData; MCIPosData];
OutboundDataAnovaGroups = [repmat({'Young'}, size(YoungData, 1), 1); ...
          repmat({'Healthy Older'}, size(HealthyControlsData, 1), 1); ...
          repmat({'MCI Unknown'}, size(MCIUnkData, 1), 1); ...
          repmat({'MCI Negative'}, size(MCINegData, 1), 1); ...
          repmat({'MCI Positive'}, size(MCIPosData, 1), 1)];

% Anova
[p, tbl, stats] = anova1(OutboundDataAnova, OutboundDataAnovaGroups, 'off');
disp("Anova on Pitch - 1 Young 2 Elderly 3 Mci unk 4 Mci neg 5 Mci pos");
tbl
% Multiple comparisons
[results, means, ~, ~] = multcompare(stats, "Alpha", 0.01, "CType","bonferroni", 'Display','off');
disp("Anova multiple comparisons results - 1 Young 2 Elderly 3 Mci unk 4 Mci neg 5 Mci pos");
results

%
currFig = figure("Position",plotInfo.FigurePosition,Visible="off");

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

hold on;

% Creating box plots
boxplots = boxplot(OutboundData,"Notch","on","symbol","","Colors",config.color_scheme_group, "Widths",0.6);
boxes = findobj(gca,'Tag','Box');
boxes = flip(boxes);
medians = findobj(gca,'Tag','Median');
medians = flip(medians);

% Adding a patch to the box plots
for j = 1:length(boxes)
    pp = patch(get(boxes(j),'XData'), get(boxes(j), 'YData'), config.color_scheme_group(j,:), 'FaceAlpha', plotInfo.PatchAlpha);
    pp.LineStyle = "none";
end
clear j;

% Changing color of the median 
for j = 1:length(medians)
    plot(medians(j).XData,medians(j).YData,Color=plotInfo.medianColor, LineStyle="-", LineWidth=plotInfo.medianWidth);
end
clear j;

% Adding the dots for each participant
for j = 1:width(OutboundData)
    sh = scatter(j*ones(height(OutboundData),1), OutboundData(:,j));
    sh.SizeData = plotInfo.MarkerSize;
    sh.MarkerEdgeColor = "none";
    sh.MarkerFaceColor = config.color_scheme_group(j,:);
    sh.MarkerFaceAlpha = plotInfo.MarkerAlpha;
end
clear j;

% Adding error bars
ax_errorBar = errorbar(xDatameans,DataMeans,DataSems);
ax_errorBar.Color = [0 0 0];
ax_errorBar.LineWidth = 3;
ax_errorBar.LineStyle = "none";

% Adding black diamond
sc_means = scatter(xDatameans,DataMeans);
sc_means.Marker = "diamond";
sc_means.SizeData = plotInfo.meanMarkerSize;
sc_means.MarkerFaceAlpha = 1;
sc_means.MarkerFaceColor = "white";
sc_means.MarkerEdgeColor = "none";

sigstaroptions.textSize      = plotInfo.sigmaStarTextSize;
sigstaroptions.lineWidth     = plotInfo.sigmaStarLineWidth;
sigstaroptions.barSeparation = plotInfo.sigmaBarSeparation;

hold off;

ylim(plotInfo.yLim);
xlim(plotInfo.xLim);

ylabel("Pitch angle (degree)");

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

exportgraphics(currFig,config.ResultFolder+"/OutboundTurnPitchAngle.png",'Resolution',300);
exportgraphics(currFig,config.ResultFolder+"/OutboundTurnPitchAngle.pdf",'Resolution',300, 'ContentType','vector');

% Final cleanup to leave workspace as the end of the Preprocessing stage.
% Remove if you want to take a look at the output data.
clearvars -except config YoungControls HealthyControls MCINeg MCIPos MCIUnk

%% ------------------------------------------------------------------
function dataout = extractPitchData(groupData)

    dataout = nan(length(groupData.Path),1);
    
    for i_participant = 1:length(groupData.Path)
        if(~isempty(groupData.Path{i_participant}))
            for i_trial = 1:length(groupData.Path{i_participant})
                path = groupData.Path{i_participant}{i_trial};
                % Including only tracking data between reaching cone 2 and
                % cone 3
                filter_idxs = path(:,1) > groupData.FlagTrigTimes{i_participant}{i_trial,2} & path(:,1) < groupData.FlagTrigTimes{i_participant}{i_trial,3};
                % Inverting the order to take into account Unity is y-up
                % coordinate system
                path = path(filter_idxs,[5 7 6]);
                curr_pitch_angle = atan2(path(:,3), sqrt(path(:,1).^2 + path(:,2).^2)) * 180 / pi;
            end   
            dataout(i_participant,1) = mean(curr_pitch_angle);
        end
    end

end

%% ------------------------------------------------------------------
function drawHeadingPoints(points, timestamps, startTime, endTime)

filter_idxs = timestamps > startTime & timestamps < endTime;
points = points(filter_idxs,:);
nPoints = height(points);
cmap    = colormap(cool);
idxs = round(linspace(1,size(cmap,1),nPoints)');
tmap    = cmap(idxs,:); 
clear cmap idxs

figure('Position', [200, 200, 1600, 800]);
subplot(1,2,1);
axis([-1, 1, -1, 1, -1, 1]);
view(45, 30);
xlabel("x");
ylabel("y");
zlabel("z");
hold on

subplot(1,2,2);
axis([0, nPoints, -180, 180]);
xlabel('Time');
ylabel('Pitch Angle (degrees)');
hold on

    % Loop through each point and plot it using drawnow
    for i = 1:nPoints
        hold on
        subplot(1,2,1);
        % Plot the current point
        plot3(points(i,1), points(i,2), points(i,3), '.', 'MarkerSize', 10, 'Color',tmap(i,:));
        % Update the plot and pause for a short amount of time
        hold on
        subplot(1,2,2);
        pitch_angle = atan2(points(i,3), sqrt(points(i,1)^2 + points(i,2)^2)) * 180 / pi;
        plot(i, pitch_angle, '.', 'MarkerSize', 10, 'Color',tmap(i,:));
        drawnow;
        pause(0.01);
    end

hold off

end