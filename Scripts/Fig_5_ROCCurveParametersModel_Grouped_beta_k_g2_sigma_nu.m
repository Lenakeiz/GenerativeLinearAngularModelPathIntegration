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
config.useTrialFilter = true;
config.ModelName        =   "beta_k_g2_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "sigma", "nu"];
config.NumParams        =   length(config.ParamName); % Set 100 here to avoid producing the model
% Run the model
VAM

%% Model run completed, preparing the data for plotting figures
config.ResultFolder = pwd + "/Output/ModelFigures/Fig5/ModelParameters/Grouped_beta_k_g2_sigma_nu";
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
close all;
plotInfo.defaultTextSize = 20;
plotInfo.defaultLineSize = 1.4;
plotInfo.titleFontSize = 16;
plotInfo.labelSize = 15;
plotInfo.axisSize = 14;
plotInfo.lineAlpha = 0.6;
plotInfo.YLabel = "True positive rate";
plotInfo.XLabel = "False positive rate";
plotInfo.Title = "Healthy controls / pooled MCI";
plotInfo.visible = "on";

disp("%%%%%%%%%%%%%%% ROC MCI pooled vs HC - model parameters estimation %%%%%%%%%%%%%%%")
rng("shuffle");
plotROCParametersCurve(HealthyControlsParameters, MCIAllParameters,'HC', 'MCI', config, plotInfo);
%
plotInfo.Title = "MCI - / MCI +";
disp("%%%%%%%%%%%%%%% ROC MCI positive vs MCI negative - model parameters estimation %%%%%%%%%%%%%%%")
rng("shuffle");
plotROCParametersCurve(MCINegParameters, MCIPosParameters,'MCIneg', 'MCIpos', config, plotInfo);

%%
function plotROCParametersCurve(params1, params2, params1groupName, params2groupName, config, plotInfo)

parametersName = [{'\beta'}, {'k'}, {'g_2'}, {'\sigma'}, {'\nu'}];
colors = config.color_scheme_npg([8 2 3 10 5 9],:);

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

AUC = plotROCCurveSVM(params1, params2, params1groupName, params2groupName, config.color_scheme_npg(4,:));

legendText{1,1} = "AUC(" + convertCharsToStrings({'\beta k g_2 \sigma \nu'}) + ") = " + num2str(round(AUC.mean,2),2);
disp("AUC std (all params): " + num2str(AUC.std));
drawnow;

AUC = plotROCCurveSVM(params1(:,[1 2 3]), params2(:,[1 2 3]), params1groupName, params2groupName, config.color_scheme_npg(8,:));

legendText{1,2} = "AUC(" + convertCharsToStrings({'\beta k g_2'}) + ") = " + num2str(round(AUC.mean,2),2);
disp("AUC std (beta k): " + num2str(AUC.std));
drawnow;

AUC = plotROCCurveSVM(params1(:,[3 4]), params2(:,[3 4]), params1groupName, params2groupName, config.color_scheme_npg(2,:));

legendText{1,3} = "AUC(" + convertCharsToStrings({'\sigma \nu'}) + ") = " + num2str(round(AUC.mean,2),2);
disp("AUC std (sigma nu): " + num2str(AUC.std));
drawnow;

hline = refline([1 0]);
hline.Color = 'k';
hline.LineWidth = 2;
hline.LineStyle = '--';

hold off

ll = legend('Location','southeast');
ll.String = legendText;
ll.FontSize = 8;

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
function AUCOut = plotROCCurveSVM(param1, param2, param1Label, param2Label, paramColor)
    
    % logistic regression with cross validation
    AUCOut=struct;
    AUCOut.mean = 0;
    AUCOut.CI   = [0 0];
    AUCOut.std  = 0;

    % Preparing the logistic regression
    allData = [param1; param2];
    allDatalogicalResponse = (1:height(param1) + height(param2))' > height(param1);
    label1     = cell(height(param1),1);
    label1(:)  = {param1Label};
    label2     = cell(height(param2),1);
    label2(:)  = {param2Label};
    allLabels  = [label1;label2];
    
    M = 1000;

    X_All = cell(1,M);
    Y_All = cell(1,M);
    AUC_All = zeros(1,M);

    holdout = 0.3;
    N = ceil(holdout*size(allData,1));
    warning('off')

    for i =1:M
        mdl  = fitcsvm(allData,allLabels,"Standardize",true,"ClassNames",[{param1Label},{param2Label}], CrossVal="on", Holdout=holdout, KernelFunction="linear");
        comp_mdl_post = fitPosterior(mdl.Trained{1},allData,allLabels);
        %predict on leave-out testing data
        [~,post_probabilities] = predict(comp_mdl_post,allData(mdl.Partition.test,:));
        [X,Y,~,AUC] = perfcurve(allDatalogicalResponse(mdl.Partition.test,:), post_probabilities(:,2), true);
        AUC_All(i) = AUC(1);

        if size(X,1)<N
            % pad 1
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

    end
    
    warning('on')
    X_mean = mean(cell2mat(X_All),2,"omitnan"); %X_std = std(New_X_All,0,2);
    Y_mean = mean(cell2mat(Y_All),2,"omitnan"); Y_std = std(cell2mat(Y_All),0,2,"omitnan")./sqrt(N);
        
    AUCOut.mean = mean(AUC_All,2,"omitnan");
    AUCOut.CI   = [0 0];
    AUCOut.std  = std(AUC_All,0,"omitnan");

    hold on
    plot(X_mean,Y_mean,"Color",paramColor,'LineWidth',2.0);
    %plot(intervals,Y_mean_manual,"Color",[paramColor 0.2],'LineWidth',2.0);
    patch([X_mean; flipud(X_mean)], [Y_mean+Y_std; flipud(Y_mean-Y_std)], paramColor, 'EdgeColor','none', 'FaceAlpha',0.2, 'HandleVisibility','off')
    hold off

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