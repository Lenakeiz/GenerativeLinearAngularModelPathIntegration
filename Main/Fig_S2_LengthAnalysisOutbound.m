%% Script to create output for Supp Fig. 2 - analysis of triangle lenghts
% Andrea Castegnaro, UCL, 2022 andrea.castegnaro@ucl.ac.uk
% Calculate length of outbound paths and visualize them in boxplots for
% different groups

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
config.ResultFolder = pwd + "/Output/S2/LengthAnalysis";
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Generating color scheme for our paper
ColorPattern;

%% Collecting information from output
YoungOutbound = extractLengthData(YoungControls);
HealthyControlsOutbound = extractLengthData(HealthyControls);
MCINegOutbound = extractLengthData(MCINeg);
MCIPosOutbound = extractLengthData(MCIPos);
MCIUnkOutbound = extractLengthData(MCIUnk);

% Separate MCI
close all;

% Parameters set for controlling visual output
plotInfo.defaultLineSize = 1.7;
plotInfo.titleFontSize = 12;
plotInfo.labelSize = 12;
plotInfo.axisSize = 10;
plotInfo.MarkerSize = 20;
plotInfo.MarkerAlpha = 0.35;
plotInfo.PatchAlpha = 0.7;
plotInfo.yLim = [4 8.5];
plotInfo.xLim = [0.5 5.5];
plotInfo.medianColor = [0.4 0.4 0.4];
plotInfo.medianWidth = 1.3;
plotInfo.meanMarkerSize = 30;
plotInfo.sigmaStarLineWidth = 2.5;
plotInfo.sigmaStarTextSize  = 20;
plotInfo.sigmaBarSeparation = 0.04;
plotInfo.FigurePosition = [200 200 280 250];


% Prepare data for plotting
max_element = max([numel(YoungOutbound),...
    numel(HealthyControlsOutbound),...
    numel(MCINegOutbound),...
    numel(MCIPosOutbound),...
    numel(MCIUnkOutbound)]);
YoungOutbound(end+1 : max_element) = nan;
HealthyControlsOutbound(end+1 : max_element) = nan;
MCIUnkOutbound(end+1 : max_element) = nan;
MCINegOutbound(end+1 : max_element) = nan;
MCIPosOutbound(end+1 : max_element) = nan;

OutboundData = [YoungOutbound HealthyControlsOutbound MCIUnkOutbound MCINegOutbound MCIPosOutbound];

DataMeans = mean(OutboundData,1,"omitnan");
DataSems = std(OutboundData,1,"omitnan")./sqrt(sum(~isnan(OutboundData),1));

xDatameans = 1:length(DataMeans);

% Calculating one-way anova to assess difference mean between groups
OutboundDataAnova = [YoungOutbound; HealthyControlsOutbound; MCIUnkOutbound; MCINegOutbound; MCIPosOutbound];
OutboundDataAnovaGroups = [repmat({'Young'}, size(YoungOutbound, 1), 1); ...
          repmat({'Healthy Older'}, size(HealthyControlsOutbound, 1), 1); ...
          repmat({'MCI Unknown'}, size(MCIUnkOutbound, 1), 1); ...
          repmat({'MCI Negative'}, size(MCINegOutbound, 1), 1); ...
          repmat({'MCI Positive'}, size(MCIPosOutbound, 1), 1)];

% Anova
[p, tbl, stats] = anova1(OutboundDataAnova, OutboundDataAnovaGroups, 'off');
disp("Anova on Length- 1 Young 2 Elderly 3 Mci unk 4 Mci neg 5 Mci pos");
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

% Adding mean value
sc_means = scatter(xDatameans,DataMeans);
sc_means.Marker = "diamond";
sc_means.SizeData = plotInfo.meanMarkerSize;
sc_means.MarkerFaceAlpha = 1;
sc_means.MarkerFaceColor = "white";
sc_means.MarkerEdgeColor = "none";

sigstaroptions.textSize      = plotInfo.sigmaStarTextSize;
sigstaroptions.lineWidth     = plotInfo.sigmaStarLineWidth;
sigstaroptions.barSeparation = plotInfo.sigmaBarSeparation;

% Manually set from looking at results table (please see above)
adjustablesigstar([1 2],0.001,0,sigstaroptions);

hold off;

ylim(plotInfo.yLim);
xlim(plotInfo.xLim);

ylabel("Length (m)");

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

exportgraphics(currFig,config.ResultFolder+"/OutboundPathLength.png",'Resolution',300);
exportgraphics(currFig,config.ResultFolder+"/OutboundPathLength.pdf",'Resolution',300, 'ContentType','vector');

%%
function dataout = extractLengthData(groupData)

    dataout = [];
    
    for i_participant = 1:length(groupData.Results.DX)
        if(~isempty(groupData.Results.DX{i_participant}))
            paths = cell2mat(groupData.Results.DX{i_participant});
            % We are interested in the oubound path only
            paths = paths([1 2],:);
            outbound_lengths = sum(paths,1);
            outbound_mean = mean(outbound_lengths,'omitnan');
            dataout = [dataout; outbound_mean];
        end
    end
end