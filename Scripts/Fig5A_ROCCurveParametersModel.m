%% Preparing the data
VAM

%% Setting colors for using in plots
ColorPattern;
resultfolder = savefolder + "PaperFigs/Fig5";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

clear resultfolder
%% Aggregating the MCI sample
AggregateMCI;

%% Creating the roc curves for MCI all vs HC

close all;

parametersName = [{'\gamma'},  {'g_3'}, {'b'}, {'\sigma'}, {'\nu'}];

% set figure info
%f = figure('visible','off','Position', [100 100 1000 500]);
f = figure('visible','off','Position', [100 100 500 500]);
%%% Font type and size setting %%%
% Using Arial as default because all journals normally require the font to
% be either Arial or Helvetica
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12)

hold on;

AUC{1} = plotROCCurve(HealthyControlsParameters, MCIAllParameters, {'\gamma g3 b \sigma \nu'}, 'AllParamsMCIvsHC','HC', 'MCI', config.color_scheme_npg(4,:), config);
legendText{1,1} = "AUC(" + convertCharsToStrings({'\gamma g3 b \sigma \nu'}) + ") = " + num2str(round(AUC{1}.Value(1),2),2);
% % Removing values that have hitted the boudary (b == 0)
% allData = [HealthyControlsParameters; MCIAllParameters];
% allDataLogicalResponse = (1:height(HealthyControlsParameters) + height(MCIAllParameters))' > height(HealthyControlsParameters);
% label1     = cell(height(HealthyControlsParameters),1);
% label1(:)  = {'HC'};
% label2     = cell(height(MCIAllParameters),1);
% label2(:)  = {'MCI'};
% allLabels  = [label1;label2];
% 
% clear label1 label2

filter = ~(HealthyControlsParameters(:,3) == 0);
HealthyControlsParametersFilter = HealthyControlsParameters(filter,:);

filter = ~(MCIAllParameters(:,3) == 0);
MCIAllParametersFilter = MCIAllParameters(filter,:);

clear filter

parametersName = [{'\gamma'},  {'g3'}, {'b'}, {'\sigma'}, {'\nu'}];
filesName = [{'GammaMCIvsHC'},  {'g3MCIvsHC'}, {'bMCIvsHC'}, {'SigmaMCIvsHC'}, {'NuMCIvsHC'}];
colors = config.color_scheme_npg([8 3 7 9 10],:);

for i = 1:length(parametersName)
    if(i == 3)
        % means b which has a lot of zeros
        AUC{i + 1}=plotROCCurve(HealthyControlsParametersFilter(:,i), MCIAllParametersFilter(:,i), parametersName{i}, filesName{i}, 'HC', 'MCI', colors(i,:), config);
    else
        AUC{i + 1}=plotROCCurve(HealthyControlsParameters(:,i), MCIAllParameters(:,i), parametersName{i}, filesName{i}, 'HC', 'MCI', colors(i,:), config);
    end
    legendText{1,i + 1} = "AUC(" + convertCharsToStrings(parametersName{i}) + ") = " + num2str(round(AUC{i + 1}.Value(1),2),2);
end

hold off;
    
% TextH = annotation('textbox', 'Position', [0.7 0.6 0.5 0.1], ...
%     'String', {"AUC(" + convertCharsToStrings({'\gamma g_3 b \sigma \nu'}) + ") = " + num2str(round(AUC{1}.Value(1),2),2)}, ...
%     'HorizontalAlignment', 'left', ...
%     'VerticalAlignment', 'top',...
%     'FitBoxToText','on', ...
%     'FontSize', 12);

ll = legend('Location','southeast');
ll.String = legendText;
ll.FontSize = 12;

ylabel('True Positive Rate');
xlabel('False Positive Rate')
title('MCI Merged / Healthy Old controls - Base Model');

%Further post-processing the figure
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...
    'LineWidth'   , .5        );

axis square;

exportgraphics(f,config.ResultFolder+"/ROC_MCIvsHC"+".png",'Resolution',300);
exportgraphics(f,config.ResultFolder+"/ROC_MCIvsHC"+".pdf",'Resolution',300, 'ContentType','vector');

clear parametersName filesName i colors f ll legendText 

%% Creating the roc curves for MCIpos vs MCIneg
close all;

parametersName = [{'\gamma'},  {'g3'}, {'b'}, {'\sigma'}, {'\nu'}];

% set figure info
%f = figure('visible','off','Position', [100 100 1000 500]);
f = figure('visible','off','Position', [100 100 500 500]);
%%% Font type and size setting %%%
% Using Arial as default because all journals normally require the font to
% be either Arial or Helvetica
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12)

hold on;

AUC{1} = plotROCCurve(MCINegParameters, MCIPosParameters, {'\gamma g3 b \sigma \nu'}, 'AllParamsMCIvsHC','MCIneg', 'MCIpos', config.color_scheme_npg(4,:), config);
legendText{1,1} = "AUC(" + convertCharsToStrings({'\gamma g3 b \sigma \nu'}) + ") = " + num2str(round(AUC{1}.Value(1),2),2);
% % Removing values that have hitted the boudary (b == 0)
% allData = [HealthyControlsParameters; MCIAllParameters];
% allDataLogicalResponse = (1:height(HealthyControlsParameters) + height(MCIAllParameters))' > height(HealthyControlsParameters);
% label1     = cell(height(HealthyControlsParameters),1);
% label1(:)  = {'HC'};
% label2     = cell(height(MCIAllParameters),1);
% label2(:)  = {'MCI'};
% allLabels  = [label1;label2];
% 
% clear label1 label2

filter = ~(HealthyControlsParameters(:,3) == 0);
HealthyControlsParametersFilter = HealthyControlsParameters(filter,:);

filter = ~(MCIAllParameters(:,3) == 0);
MCIAllParametersFilter = MCIAllParameters(filter,:);

clear filter

parametersName = [{'\gamma'},  {'g3'}, {'b'}, {'\sigma'}, {'\nu'}];
filesName = [{'GammaMCIvsHC'},  {'g3MCIvsHC'}, {'bMCIvsHC'}, {'SigmaMCIvsHC'}, {'NuMCIvsHC'}];
colors = config.color_scheme_npg([8 3 7 9 10],:);

for i = 1:length(parametersName)
    if(i == 3)
        % means b which has a lot of zeros
        AUC{i + 1}=plotROCCurve(MCINegParameters(:,i), MCIPosParameters(:,i), parametersName{i}, filesName{i}, 'MCIneg', 'MCIpos', colors(i,:), config);
    else
        AUC{i + 1}=plotROCCurve(MCINegParameters(:,i), MCIPosParameters(:,i), parametersName{i}, filesName{i}, 'MCIneg', 'MCIpos', colors(i,:), config);
    end
    legendText{1,i + 1} = "AUC(" + convertCharsToStrings(parametersName{i}) + ") = " + num2str(round(AUC{i + 1}.Value(1),2),2);
end

hold off;
    
% TextH = annotation('textbox', 'Position', [0.7 0.6 0.5 0.1], ...
%     'String', {"AUC(" + convertCharsToStrings({'\gamma g_3 b \sigma \nu'}) + ") = " + num2str(round(AUC{1}.Value(1),2),2)}, ...
%     'HorizontalAlignment', 'left', ...
%     'VerticalAlignment', 'top',...
%     'FitBoxToText','on', ...
%     'FontSize', 12);

ll = legend('Location','southeast');
ll.String = legendText;
ll.FontSize = 12;

ylabel('True Positive Rate');
xlabel('False Positive Rate')
title('MCI negative / MCI positive - Base Model');

%Further post-processing the figure
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...
    'LineWidth'   , .5        );

axis square;

exportgraphics(f,config.ResultFolder+"/ROC_MCInegvsMCIPos"+".png",'Resolution',300);
exportgraphics(f,config.ResultFolder+"/ROC_MCInegvsMCIPos"+".pdf",'Resolution',300, 'ContentType','vector');

clear parametersName filesName i colors f ll legendText 

%%
function AUC = plotROCCurve(param1, param2, paramName, filename, param1Label, param2Label, paramColor, config)
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
    [X,Y,~,AUC.Value] = perfcurve(allLabels, allDataScores, param2Label, NBoot=10);
    %[~,AUC.CI] = auc([allLabels allDataScores],0.05,'boot',10000,'type','bca');

    plot(X(:,1),Y(:,1),...
        "Color",paramColor,...
        'LineWidth',1.0...
        );
end