
%% Model run completed, preparing the data for plotting figures
config.ResultFolder = pwd + "/Output/ModelFigures/ROC_FH_M1";
% Create storing folder for trajectory if not exist
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

%% Genarating color scheme
ColorPattern;
%% Getting Information from results:
FHNegParameters          = averageAcrossConditions(FamilyHistNeg.Results.estimatedParams);
FHPosParameters          = averageAcrossConditions(FamilyHistPos.Results.estimatedParams); 

%% Plotting roc curve HC vs pooled MCI and MCI negative vs MCI positive
% Plotting variables
plotInfo.defaultTextSize = 20;
plotInfo.defaultLineSize = 1.4;
plotInfo.titleFontSize = 16;
plotInfo.labelSize = 15;
plotInfo.axisSize = 14;
plotInfo.lineAlpha = 0.6;
plotInfo.YLabel = "True positive rate";
plotInfo.XLabel = "False positive rate";
plotInfo.visible = "on";

plotInfo.Title = "FH negative / FH positive";

disp("%%%%%%%%%%%%%%% ROC FH positive vs FH negative - model parameters estimation %%%%%%%%%%%%%%%")
plotROCParametersCurve(FHNegParameters, FHPosParameters,'MCIneg', 'MCIpos', config, plotInfo);

%%
function plotROCParametersCurve(params1, params2, params1groupName, params2groupName, config, plotInfo)

parametersName = [{'\sigma'}, {'\nu'}];
colors = config.color_scheme_npg([8 3 7 9 10],:);

% set figure info
%f = figure('visible','off','Position', [100 100 1000 500]);
f = figure('visible', plotInfo.visible, 'Position', [0 0 500 400]);
%%% Font type and size setting %%%
% Using Arial as default because all journals normally require the font to
% be either Arial or Helvetica
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',plotInfo.axisSize)
set(0,'DefaultTextFontSize',plotInfo.axisSize)

hold on;

AUC{1} = plotsingleROCCurve(params1, params2, params1groupName, params2groupName, config.color_scheme_npg(4,:), plotInfo);
legendText{1,1} = "AUC(" + convertCharsToStrings({'\sigma \nu'}) + ") = " + num2str(round(AUC{1}.Value(1),2),2);

disp("AUC(" + convertCharsToStrings({'\sigma \nu'}) + ") CI = [" ...
    + num2str(round(AUC{1}.Value(2),2),2) +  ", " + num2str(round(AUC{1}.Value(3),2),2) + "]")
% filter = ~(HealthyControlsParameters(:,3) == 0);
% HealthyControlsParametersFilter = HealthyControlsParameters(filter,:);
% 
% filter = ~(MCIAllParameters(:,3) == 0);
% MCIAllParametersFilter = MCIAllParameters(filter,:);

% clear filter

for i = 1:length(parametersName)
    AUC{i + 1} = plotsingleROCCurve(params1(:,i), params2(:,i), params1groupName, params2groupName, colors(i,:), plotInfo);
    legendText{1,i + 1} = "AUC(" + convertCharsToStrings(parametersName{i}) + ") = " + num2str(round(AUC{i + 1}.Value(1),2),2);
    disp("AUC(" + convertCharsToStrings(parametersName{i}) + ") CI = [" ...
        + num2str(round(AUC{i + 1}.Value(2),2),2) +  ", " + num2str(round(AUC{i + 1}.Value(3),2),2) + "]")
end

hold off;

ll = legend('Location','southeast');
ll.String = legendText;
ll.FontSize = 12;

ylabel('True Positive Rate');
xlabel('False Positive Rate')
t = title(plotInfo.Title);

t.FontSize = plotInfo.titleFontSize;

%Further post-processing the figure
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...
    'LineWidth'   , .5        );

axis square;

ax = gca;
ax.LineWidth = plotInfo.defaultLineSize;
ax.XLabel.FontSize = plotInfo.labelSize;
ax.YLabel.FontSize = plotInfo.labelSize;
ax.XAxis.FontSize = plotInfo.axisSize;
ax.YAxis.FontSize = plotInfo.axisSize;
ax.YTick = 0:0.2:1.0;

exportName = [params1groupName 'vs' params2groupName];

exportgraphics(f,config.ResultFolder+"/"+convertCharsToStrings(exportName)+".png",'Resolution',300);
exportgraphics(f,config.ResultFolder+"/"+convertCharsToStrings(exportName)+".pdf",'Resolution',300, 'ContentType','vector');

clear parametersName filesName i colors f ll legendText
end

%%
function AUC = plotsingleROCCurve(param1, param2, param1Label, param2Label, paramColor, plotInfo)
    % Preparing the logistic regression
    allData = [param1; param2];
    allDatalogicalResponse = (1:height(param1) + height(param2))' > height(param1);
    label1     = cell(height(param1),1);
    label1(:)  = {param1Label};
    label2     = cell(height(param2),1);
    label2(:)  = {param2Label};
    allLabels  = [label1;label2];

    clear hcLabel mciLabel
    
    % Fitting the logistic regression
    mdl = fitglm(allData,allDatalogicalResponse,'Distribution', 'binomial','Link','logit');
    
    allDataScores = mdl.Fitted.Probability;
    [X,Y,~,AUC.Value] = perfcurve(allLabels, allDataScores, param2Label, NBoot=100);

    lineplot = plot(X(:,1),Y(:,1),...
        "Color",paramColor,...
        'LineWidth',2.0...
        );

    lineplot.Color = [lineplot.Color plotInfo.lineAlpha];

end

%% get the model parameters, average across the conditions
% remove nans if the row after mergin still contains nans
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