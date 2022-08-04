%% Preparing the data
VAM_PrepareBaseConfig

%% Preprocessing the data
VAM_PreprocessData

%% Setting the model we are interested in
% Eventually modify config paramteters we are interested in. For example
% for this graph we are not interested in running the model with splitted
% conditions so we will set the relative config to false 
% force it to not run
config.useTrialFilter = true;
config.ModelName        =   "beta_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName); % Set 100 here to avoid producing the model
% Run the model
VAM

%% Model run completed, preparing the data for plotting figures
config.ResultFolder = pwd + "/Output/PaperFigs/Fig5A";
% Create storing folder for trajectory if not exist
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

%% Genarating color scheme
ColorPattern;
%% Getting Information from results:
YoungControlsParameters   = averageAcrossConditions(YoungControls.Results.estimatedParams);
HealthyControlsParameters = averageAcrossConditions(HealthyControls.Results.estimatedParams);
MCIUnkParameters          = averageAcrossConditions(MCIUnk.Results.estimatedParams);
MCINegParameters          = averageAcrossConditions(MCINeg.Results.estimatedParams);
MCIPosParameters          = averageAcrossConditions(MCIPos.Results.estimatedParams);
MCIAllParameters          = [MCIUnkParameters; MCINegParameters; MCIPosParameters]; 

%% Plotting roc curve HC vs pooled MCI and MCI negative vs MCI positive
% Plotting variables
plotInfo.defaultTextSize = 20;
plotInfo.defaultLineSize = 1.3;
plotInfo.titleFontSize = 16;
plotInfo.labelSize = 15;
plotInfo.axisSize = 14;
plotInfo.YLabel = "True positive rate";
plotInfo.XLabel = "False positive rate";
plotInfo.Title = "Healthy controls / pooled MCI";
plotInfo.visible = "off";

disp("%%%%%%%%%%%%%%% ROC MCI pooled vs HC - model parameters estimation %%%%%%%%%%%%%%%")
plotROCParametersCurve(HealthyControlsParameters, MCIAllParameters,'HC', 'MCI', config, plotInfo);

plotInfo.Title = "MCI negative / MCI positive";

disp("%%%%%%%%%%%%%%% ROC MCI positive vs MCI negative - model parameters estimation %%%%%%%%%%%%%%%")
plotROCParametersCurve(MCINegParameters, MCIPosParameters,'MCIneg', 'MCIpos', config, plotInfo);

%%
function plotROCParametersCurve(params1, params2, params1groupName, params2groupName, config, plotInfo)

parametersName = [{'\beta'},  {'g_2'}, {'g_3'}, {'\sigma'}, {'\nu'}];
colors = config.color_scheme_npg([8 3 7 9 10],:);

% set figure info
%f = figure('visible','off','Position', [100 100 1000 500]);
f = figure('visible',plotInfo.visible,'Position', [100 100 600 600]);
%%% Font type and size setting %%%
% Using Arial as default because all journals normally require the font to
% be either Arial or Helvetica
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12)

hold on;

AUC{1} = plotsingleROCCurve(params1, params2, params1groupName, params2groupName, config.color_scheme_npg(4,:));
legendText{1,1} = "AUC(" + convertCharsToStrings({'\beta g_2 g_3 \sigma \nu'}) + ") = " + num2str(round(AUC{1}.Value(1),2),2);

disp("AUC(" + convertCharsToStrings({'\beta g_2 g_3 \sigma \nu'}) + ") CI = [" ...
    + num2str(round(AUC{1}.Value(2),2),2) +  ", " + num2str(round(AUC{1}.Value(3),2),2) + "]")
% filter = ~(HealthyControlsParameters(:,3) == 0);
% HealthyControlsParametersFilter = HealthyControlsParameters(filter,:);
% 
% filter = ~(MCIAllParameters(:,3) == 0);
% MCIAllParametersFilter = MCIAllParameters(filter,:);

% clear filter

for i = 1:length(parametersName)
    AUC{i + 1} = plotsingleROCCurve(params1(:,i), params2(:,i), params1groupName, params2groupName, colors(i,:));
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
function AUC = plotsingleROCCurve(param1, param2, param1Label, param2Label, paramColor)
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

    plot(X(:,1),Y(:,1),...
        "Color",paramColor,...
        'LineWidth',1.0...
        );
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