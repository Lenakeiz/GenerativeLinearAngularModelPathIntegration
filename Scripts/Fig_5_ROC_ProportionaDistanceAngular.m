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
config.NumParams        =   100; %length(config.ParamName); % Set 100 here to avoid producing the model
% Run the model
VAM

%% Model run completed, preparing the data for plotting figures
config.ResultFolder = pwd + "/Output/ModelFigures/Fig5/ProportionalError";
% Create storing folder for trajectory if not exist
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

%% Genarating color scheme
ColorPattern;

%%
HealthyControlsDistanceError = averageAcrossConditions(HealthyControls.Results.PropDistErr);
MCIUnkDistanceError          = averageAcrossConditions(MCIUnk.Results.PropDistErr);
MCINegDistanceError          = averageAcrossConditions(MCINeg.Results.PropDistErr);
MCIPosDistanceError          = averageAcrossConditions(MCIPos.Results.PropDistErr);
MCIAllDistanceError          = [MCIUnkDistanceError; MCINegDistanceError; MCIPosDistanceError]; 

HealthyControlsAngErr        = averageAcrossConditions(HealthyControls.Results.PropAngErr);
MCIUnkAngErr                 = averageAcrossConditions(MCIUnk.Results.PropAngErr);
MCINegAngErr                 = averageAcrossConditions(MCINeg.Results.PropAngErr);
MCIPosAngErr                 = averageAcrossConditions(MCIPos.Results.PropAngErr);
MCIAllAngErr                 = [MCIUnkAngErr; MCINegAngErr; MCIPosAngErr];

%
allParamsHC        = [HealthyControlsDistanceError HealthyControlsAngErr];
allParamsPooledMCI = [MCIAllDistanceError MCIAllAngErr];
allParamsMCIPos    = [MCIPosDistanceError MCIPosAngErr];
allParamsMCINeg    = [MCINegDistanceError MCINegAngErr];

%% Plotting roc curve HC vs pooled MCI and MCI negative vs MCI positive
% Plotting variables
close all;

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
%%
plotInfo.Title = "MCI negative / MCI positive";
disp("%%%%%%%%%%%%%%% ROC MCI positive vs MCI negative - behavioural data %%%%%%%%%%%%%%%");
rng("default");
generateROCCurve(allParamsMCINeg, allParamsMCIPos,'MCI-', 'MCI+', parametersName, config, plotInfo);

%%
function generateROCCurve(params1, params2, params1groupName, params2groupName, parametersName, config, plotInfo)

colors = config.color_scheme_npg([8 3 6 9 2 4],:);
% set figure info
f = figure('visible', plotInfo.visible, 'Position', [200 200 250 200]);
%%% Font type and size setting %%%
% Using Arial as default because all journals normally require the font to
% be either Arial or Helvetica
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',plotInfo.axisSize)
set(0,'DefaultTextFontSize',plotInfo.axisSize)

hold on;

for i = 1:length(parametersName)
    %AUC = plotsingleROCCurve(params1(:,i), params2(:,i), params1groupName, params2groupName, colors(i,:));
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
%t = title(plotInfo.Title);
%t.FontSize = plotInfo.titleFontSize;

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
function AUCOut = plotsingleROCCurveSVM(param1, param2, param1Label, param2Label, paramColor)
    % logistic regression with cross validation
    AUCOut=struct;
    % Preparing the logistic regression
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
    %cross validation for M times
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
        %predict on leave-out testing data
        [predicted_labels,post_probabilities] = predict(comp_mdl_post,Xtest);
        [X,Y,~,AUC] = perfcurve(Ytest, post_probabilities(:,2), {param2Label});
        N = length(Ytest);

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
    %plot(intervals,Y_mean_manual,"Color",[paramColor 0.2],'LineWidth',2.0);
    patch([X_mean; flipud(X_mean)], [Y_mean+Y_std; flipud(Y_mean-Y_std)], paramColor, 'EdgeColor','none', 'FaceAlpha',0.2, 'HandleVisibility','off')
    hold off
end

function x= adjust_unique_points(Xroc)
    x = zeros(size(Xroc,1), size(Xroc,2));
    aux= 0.000001;
    for i=1:size(Xroc,1)
        if i~=1
            x(i,:)= Xroc(i,:)+aux;
            aux=aux+0.0001;
        end        
    end
end

%%
function AUCOut = plotsingleROCCurve(param1, param2, param1Label, param2Label, paramColor)
    % logistic regression with cross validation
    
    AUCOut=struct;

    % Preparing the logistic regression
    allData = [param1; param2];
    allDatalogicalResponse = (1:height(param1) + height(param2))' > height(param1);
    label1     = cell(height(param1),1);
    label1(:)  = {param1Label};
    label2     = cell(height(param2),1);
    label2(:)  = {param2Label};
    allLabels  = [label1;label2];

    M=1000;
    propLeaveOut = 0.3;
    N = ceil(propLeaveOut*size(allData,1));
    X_All = cell(1,M);
    Y_All = cell(1,M);

    X_All_low = cell(1,M);
    Y_All_low = cell(1,M);

    X_All_high = cell(1,M);
    Y_All_high = cell(1,M);

    AUC_All = zeros(1,M);
    AUC_CI = zeros(M,2);

    warning('off',"all");
    %cross validation for M times
    for i=1:M
        cv1 = cvpartition(size(param1,1),'HoldOut',propLeaveOut); %leave 30% datat out for the first group
        cv2 = cvpartition(size(param2,1),'HoldOut',propLeaveOut); % leave 30% data out for the second group
        train_idx = [cv1.training;cv2.training];
        trainData = allData(train_idx,:);
        trainLabel = allDatalogicalResponse(train_idx,:);
        
        test_idx = [cv1.test;cv2.test];
        testData = allData(test_idx,:);

        %testLabel = allDatalogicalResponse(cv.test,:);

        % Fitting the logistic regression
        mdl = fitglm(trainData,trainLabel,'Distribution', 'binomial','Link','logit');
        
        %predict on leave-out testing data
        testScore = predict(mdl, testData);

        testLabels = allLabels(test_idx,:);

        [X,Y,~,AUC] = perfcurve(testLabels, testScore, param2Label, NBoot=100);

        if isnan(X(2:3))
            AUC_All(i) = nan;
            AUC_CI(i,:)  = nan(1,2);
            continue
        end

        if length(X)<N
            %pad 1 after the vector
            %not sure why length of elements in X_All and Y_All are not consistent
            %so far. May be due to ROC calculation somewhere inside the code
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

        X_All_low = X(:,2);
        Y_All_low = Y(:,2);

        X_All_high = X(:,3);
        Y_All_high = Y(:,3);

        AUC_All(i) = AUC(1);
        AUC_CI(i,:)  = AUC(2:3);
    end

    X_mean = mean(cell2mat(X_All),2);
    Y_mean = mean(cell2mat(Y_All),2); Y_d = std(cell2mat(Y_All),0,2);%./sqrt(M);

    AUCOut.mean = mean(AUC_All,2,"omitnan");
    AUCOut.CI   = mean(AUC_CI,1,"omitnan");
    warning('on',"all");
    hold on;
    plot(X_mean,Y_mean, "Color",paramColor,'LineWidth',2.0);
    patch([X_mean; flipud(X_mean)], [Y_mean+Y_d; flipud(Y_mean-Y_d)], paramColor, 'EdgeColor','none', 'FaceAlpha',0.2, 'HandleVisibility','off')
    hold off;

end


%% get the model parameters, average across the conditions
% remove nans if the row after mergin still contains nans
function dataout = averageAcrossConditions(data)

    dataout = [];
    pSize   = length(data{1});
    condSize = width(data);

    for i = 1:pSize
        tempP = [];
        for j = 1:condSize % Trial conditions
            tempP = [tempP cell2mat(data{j}{i})];
        end
        tempavg = mean(tempP,"omitnan");
        dataout = [dataout;tempavg];
    end

    dataout = removeNanRows(dataout);

end