%% Preparing the data
VAM_PrepareBaseConfig

%% Preprocessing the data
VAM_PreprocessData
 
%% Setting the model we are interested in
% Eventually modify config paramteters we are interested in. For example
% for this graph we are not interested in running the model with splitted
% conditions so we will set the relative config to false 
% force it to not run
rng("default");
config.useTrialFilter = false;
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   100; % Set 100 here to avoid producing the model
% Run the model
VAM

%% Model run completed, preparing the data for plotting figures
config.ResultFolder = pwd + "/Output/Supplementary/LengthAnalysis";
% Create storing folder for trajectory if not exist
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

%% Genarating color scheme
ColorPattern;

%% Getting speed information
% Averaging across segments L1 and L2 of the outbound path
YoungOutbound = extractLengthData(YoungControls);
HealthyControlsOutbound = extractLengthData(HealthyControls);
MCINegOutbound = extractLengthData(MCINeg);
MCIPosOutbound = extractLengthData(MCIPos);
MCIUnkOutbound = extractLengthData(MCIUnk);
%%
close all;

plotInfo.defaultLineSize = 1.7;
plotInfo.titleFontSize = 12;
plotInfo.labelSize = 12;
plotInfo.axisSize = 10;
plotInfo.MarkerSize = 70;
plotInfo.MarkerAlpha = 0.3;
plotInfo.PatchAlpha = 0.7;
plotInfo.yLim = [4 8];
plotInfo.xLim = [0.5 5.5];
plotInfo.medianColor = [0.4 0.4 0.4];
plotInfo.medianWidth = 1.3;
plotInfo.meanMarkerSize = 90;
plotInfo.FigurePosition = [200 200 600 300];

% Generate example data
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

OutboundDurationData = [YoungOutbound HealthyControlsOutbound MCIUnkOutbound MCINegOutbound MCIPosOutbound];

DataMeans = mean(OutboundDurationData,1,"omitnan");
DataSems = std(OutboundDurationData,1,"omitnan")./sqrt(sum(~isnan(OutboundDurationData),1));

xDatameans = 1:length(DataMeans);

currFig = figure("Position",plotInfo.FigurePosition);

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

hold on;

% Create box plot
boxplots = boxplot(OutboundDurationData,"Notch","on","symbol","","Colors",config.color_scheme_group, "Widths",0.6);
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

for j = 1:width(OutboundDurationData)
    sh = scatter(j*ones(height(OutboundDurationData),1), OutboundDurationData(:,j));
    sh.SizeData = plotInfo.MarkerSize;
    sh.MarkerEdgeColor = "none";
    sh.MarkerFaceColor = config.color_scheme_group(j,:);
    sh.MarkerFaceAlpha = plotInfo.MarkerAlpha;
end
clear j;

hold off;

ylim(plotInfo.yLim);
xlim(plotInfo.xLim);

ylabel("Outbound path length (m)");

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

%
OutboundDurationDataAnova = [YoungOutbound; HealthyControlsOutbound; MCIUnkOutbound; MCINegOutbound; MCIPosOutbound];
OutboundDurationDataAnovaGroups = [repmat({'Young'}, size(YoungOutbound, 1), 1); ...
          repmat({'Healthy Older'}, size(HealthyControlsOutbound, 1), 1); ...
          repmat({'MCI Unknown'}, size(MCIUnkOutbound, 1), 1); ...
          repmat({'MCI Negative'}, size(MCINegOutbound, 1), 1); ...
          repmat({'MCI Positive'}, size(MCIPosOutbound, 1), 1)];

[p, tbl, stats] = anova1(OutboundDurationDataAnova, OutboundDurationDataAnovaGroups, 'off');
[results, means, ~, ~] = multcompare(stats, 'alpha', 0.05,'Display','off');

exportgraphics(currFig,config.ResultFolder+"/SupplementaryLength.png",'Resolution',300);
exportgraphics(currFig,config.ResultFolder+"/SupplementaryLength.pdf",'Resolution',300, 'ContentType','vector');


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

%% Extracting speed information for each of the walked segment
%
function dataout = averageAcrossConditions(data)

    dataout = [];
    pSize = length(data{1});
    paramsSize = width(data{1});

    for i = 1:pSize
        tempP = [];
        for j = 1:paramsSize
            tempP = [tempP mean([data{1}(i,j) data{2}(i,j) data{3}(i,j)],"omitnan")];
        end
        dataout = [dataout;tempP];
    end

    dataout = removeNanRows(dataout);

end