%% Script to create output for Fig. 5 - receiver operating characteristic curve for proportional errors 
% Andrea Castegnaro, UCL, 2022 andrea.castegnaro@ucl.ac.uk
% Using support vector machine (SVM) to classify two groups based on their
% fitted parameters using an holdout strategy (60% training set 40% testing
% set). Receiver operating characteristic curves are calculated on
% posterior probabilities of the model applied to the testing set. Cross
% validation has been applied using 1000 repetitions.

% Preparing the data
GLAMPI_PrepareBaseConfig

% Preprocessing the data
GLAMPI_PreprocessData

% Model fitting
config.useTrialFilter = true;
config.ModelName        =   "beta_k_g2_g3_m3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "m3", "sigma", "nu"];
% We are still fitting the model, as to allow fair comparisons we remove
% participants that do not have fitted params.
config.NumParams        =   length(config.ParamName);

GLAMPI

% Preparing output
config.ResultFolder = pwd + "/Output/Fig5/ProportionalError";
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Generating color scheme for our paper
ColorPattern;

%% Collecting information from output
% To allow fair comparisons wer are going to remove participants for which
% our model did not fit. See Methods for details on participants were excluded from the analysis.

HealthyControlsDistanceError = averageAcrossConditions(HealthyControls.Results.PropDistErr, HealthyControls.Results.estimatedParams);
MCIUnkDistanceError          = averageAcrossConditions(MCIUnk.Results.PropDistErr, MCIUnk.Results.estimatedParams);
MCINegDistanceError          = averageAcrossConditions(MCINeg.Results.PropDistErr, MCINeg.Results.estimatedParams);
MCIPosDistanceError          = averageAcrossConditions(MCIPos.Results.PropDistErr, MCIPos.Results.estimatedParams);

MCIAllDistanceError          = [MCIUnkDistanceError; MCINegDistanceError; MCIPosDistanceError];

HealthyControlsAngErr        = averageAcrossConditions(HealthyControls.Results.PropAngErr, HealthyControls.Results.estimatedParams);
MCIUnkAngErr                 = averageAcrossConditions(MCIUnk.Results.PropAngErr, MCIUnk.Results.estimatedParams);
MCINegAngErr                 = averageAcrossConditions(MCINeg.Results.PropAngErr, MCINeg.Results.estimatedParams);
MCIPosAngErr                 = averageAcrossConditions(MCIPos.Results.PropAngErr, MCIPos.Results.estimatedParams);
MCIAllAngErr                 = [MCIUnkAngErr; MCINegAngErr; MCIPosAngErr];

allParamsHC        = [HealthyControlsDistanceError HealthyControlsAngErr];
allParamsPooledMCI = [MCIAllDistanceError MCIAllAngErr];
allParamsMCIPos    = [MCIPosDistanceError MCIPosAngErr];
allParamsMCINeg    = [MCINegDistanceError MCINegAngErr];

%% Fitting SVM and plotting roc curves
% HC vs pooled MCI 
% Plotting variables
close all;

% Parameters set for controlling visual output
plotInfo.defaultTextSize = 20;
plotInfo.defaultLineSize = 1.7;
plotInfo.titleFontSize = 12;
plotInfo.labelSize = 12;
plotInfo.axisSize = 10;
plotInfo.lineAlpha = 0.6;
plotInfo.YLabel = "True positive rate";
plotInfo.XLabel = "False positive rate";
plotInfo.Title = "Healthy elder / MCI";
plotInfo.visible = "on";
parametersName = ["Linear", "Angular"];

disp("%%%%%%%%%%%%%%% ROC Curve pooled MCI vs HC - behavioural data %%%%%%%%%%%%%%%")
rng("default");
generateROCCurve(allParamsHC, allParamsPooledMCI,'HC', 'MCI', parametersName, config, plotInfo);

% MCI negative vs MCI positive
plotInfo.Title = "MCI negative / MCI positive";
disp("%%%%%%%%%%%%%%% ROC MCI positive vs MCI negative - behavioural data %%%%%%%%%%%%%%%");
rng("default");
generateROCCurve(allParamsMCINeg, allParamsMCIPos,'MCI-', 'MCI+', parametersName, config, plotInfo);

disp("%%%%%%%%%%%%%%% SVM fitting completed %%%%%%%%%%%%%%%");

% Final cleanup to leave workspace as the end of the Preprocessing stage.
% Remove if you want to take a look at the output data.
%clearvars -except config YoungControls HealthyControls MCINeg MCIPos MCIUnk

%%
function generateROCCurve(params1, params2, params1groupName, params2groupName, parametersName, config, plotInfo)

colors = config.color_scheme_npg([8 3 6 9 2 4],:);

f = figure('visible', plotInfo.visible, 'Position', [200 200 250 200]);

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',plotInfo.axisSize)
set(0,'DefaultTextFontSize',plotInfo.axisSize)

hold on;

for i = 1:length(parametersName)
    AUC = plotsingleROCCurveSVM(params1(:,i), params2(:,i), params1groupName, params2groupName, colors(i,:));
    legendText{1,i} = "AUC(" + parametersName(i) + ") = "+num2str(round(AUC.mean,2));
    disp("AUC std: " + num2str(AUC.std));
    drawnow;
end

hline = refline([1 0]);
hline.Color = 'k';
hline.LineWidth = 2;
hline.LineStyle = '--';

hold off;

legend('Location','southeast')
ll = legend('Location','southeast');
ll.String = legendText;
ll.FontSize = 8;
ll.Box = "off";

ylabel('True Positive Rate');
xlabel('False Positive Rate')

%% Figure post-processing
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...
    'XLim'        , [0,1],...
    'YLim'        , [0,1],...
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
function AUCOut = plotsingleROCCurveSVM(param1, param2, param1Label, param2Label, paramColor)
    
    AUCOut=struct;
    
    % Preparing data for the svm input
    allData = [param1; param2];
    allDatalogicalResponse = (1:height(param1) + height(param2))' > height(param1);
    label1     = cell(height(param1),1);
    label1(:)  = {param1Label};
    label2     = cell(height(param2),1);
    label2(:)  = {param2Label};
    allLabels  = [label1;label2];

    M=1000;
    holdout = 0.4;
    X_All = cell(1,M);
    Y_All = cell(1,M);
    AUC_All = zeros(1,M);

    warning('off',"all");

    % Cross-validation
    for i = 1:M

        mdl  = fitcsvm(allData,allLabels,ClassNames=[{param1Label},{param2Label}],Standardize=false, CrossVal="on", Holdout=holdout, KernelFunction="linear");
        trainingIDs = training(mdl.Partition);
        Xtraining = allData(trainingIDs,:);
        Ytraining = allLabels(trainingIDs);
        testingIDs = test(mdl.Partition);
        Xtest = allData(testingIDs,:);
        Ytest = allLabels(testingIDs);

        clear trainingIDs testingIDs
        compactModel = mdl.Trained{1};

        comp_mdl_post = fitPosterior(compactModel,Xtraining,Ytraining);
        
        [predicted_labels,post_probabilities] = predict(comp_mdl_post,Xtest);
        [X,Y,~,AUC] = perfcurve(Ytest, post_probabilities(:,2), {param2Label});
        N = length(Ytest);

        if size(X,1)<N
            % pad 1 to make all of the repetition of the same length
            X = [X;ones(N-size(X,1),size(X,2))];
            Y = [Y;ones(N-size(Y,1),size(Y,2))];
        else
            X = X(1:N,:);
            X(end,:) = ones(1,size(X,2));

            Y = Y(1:N,:);
            Y(end,:) = ones(1,size(Y,2));
        end

        X_All{:,i} = X(:,1);
        Y_All{:,i} = Y(:,1);

        AUC_All(i) = AUC(1);

    end

    X_mean = mean(cell2mat(X_All),2,"omitnan");
    Y_mean = mean(cell2mat(Y_All),2,"omitnan"); 
    Y_std = std(cell2mat(Y_All),0,2,"omitnan")./sqrt(N);

    AUCOut.mean = mean(AUC_All,2,"omitnan");
    AUCOut.CI   = [0 0];
    AUCOut.std  = std(AUC_All,0,"omitnan");
    
    warning('on',"all");
    hold on
    plot(X_mean,Y_mean,"Color",paramColor,'LineWidth',2.0);
    patch([X_mean; flipud(X_mean)], [Y_mean+Y_std; flipud(Y_mean-Y_std)], paramColor, 'EdgeColor','none', 'FaceAlpha',0.2, 'HandleVisibility','off')
    hold off
end

%% 
function dataout = averageAcrossConditions(data, glampi_data)
    % average across the conditions, assuming a specific structure to this
    % data
    dataout = [];
    pSize   = length(data{1});
    paramsSize = width(data);

    for i = 1:pSize
        tempP = [];
        for j = 1:paramsSize
            tempP = [tempP cell2mat(data{j}{i})];
        end
        tempavg = mean(tempP,"omitnan");
        dataout = [dataout;tempavg];
    end

    % remove nans if the row after merging data
    % dataout = removeNanRows(dataout);
    % remove participants for which we were unable to fit the model
    excludedParticipants = removeExcludedParticipants(glampi_data);
    dataout = dataout(excludedParticipants,:);
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