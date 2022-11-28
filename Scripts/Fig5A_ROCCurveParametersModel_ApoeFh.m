
%% Model run completed, preparing the data for plotting figures
config.ResultFolder = pwd + "/Output/ModelFigures/ROC_ApoeFh";
% Create storing folder for trajectory if not exist
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

%% Genarating color scheme
ColorPattern;
%% Getting Information from results:
DoubleNegParameters          = averageAcrossConditions(DoubleNeg.Results.estimatedParams);
DoublePosParameters          = averageAcrossConditions(DoublePos.Results.estimatedParams); 

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

plotInfo.Title = "Apoe+FH Neg / Double positive";

disp("%%%%%%%%%%%%%%% ROC FH positive vs FH negative - model parameters estimation %%%%%%%%%%%%%%%")
plotROCParametersCurve(DoubleNegParameters, DoublePosParameters,'Apoe+FH Neg', 'Apoe+FH Pos', config, plotInfo);

%%
function plotROCParametersCurve(params1, params2, params1groupName, params2groupName, config, plotInfo)

parametersName = [{'\beta'},  {'g_2'}, {'g_3'}, {'\sigma'}, {'\nu'}];
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

AUC = plotsingleROCCurve(params1, params2, params1groupName, params2groupName, config.color_scheme_npg(4,:));
legendText{1,1} = "AUC(" + convertCharsToStrings({'\beta g_2 g_3 \sigma \nu'}) + ") = " + num2str(round(AUC,2),2);

for i = 1:length(parametersName)
    AUC = plotsingleROCCurve(params1(:,i), params2(:,i), params1groupName, params2groupName, colors(i,:));
    legendText{1,i+1} = "AUC(" + parametersName(i) + ") = "+num2str(round(AUC,2));
end

hline = refline([1 0]);
hline.Color = 'k';
hline.LineWidth = 2;
hline.LineStyle = '--';

hold off;

ll = legend('Location','southeast');
ll.String = legendText;
ll.FontSize = 10;

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
function AUC_mean = plotsingleROCCurve(param1, param2, param1Label, param2Label, paramColor)
    % logistic regression with cross validation

    % Preparing the logistic regression
    allData = [param1; param2];
    allDatalogicalResponse = (1:height(param1) + height(param2))' > height(param1);
    label1     = cell(height(param1),1);
    label1(:)  = {param1Label};
    label2     = cell(height(param2),1);
    label2(:)  = {param2Label};
    allLabels  = [label1;label2];

    M=500;
    N = ceil(0.3*size(allData,1));
    X_All = cell(1,M);
    Y_All = cell(1,M);
    AUC_All = zeros(1,M);

    %cross validation for M times
    for i=1:M
        cv1 = cvpartition(size(param1,1),'HoldOut',0.3); %leave 30% datat out for the first group
        cv2 = cvpartition(size(param2,1),'HoldOut',0.3); % leave 30% data out for the second group
        train_idx = [cv1.training;cv2.training];
        trainData = allData(train_idx,:);
        trainLabel = allDatalogicalResponse(train_idx,:);
        
        test_idx = [cv1.test;cv2.test];
        testData = allData(test_idx,:);

        %testLabel = allDatalogicalResponse(cv.test,:);

        % Fitting the logistic regression
        %mdl = fitglm(trainData,trainLabel,'Distribution', 'binomial','Link','logit');
        mdl = fitglm(trainData,trainLabel);

        %predict on leave-out testing data
        testScore = predict(mdl, testData);

        testLabels = allLabels(test_idx,:);

        [X,Y,~,AUC] = perfcurve(testLabels, testScore, param2Label);
        
        if length(X)<N
            %pad 1 after the vector
            %not sure why length of elements in X_All and Y_All are not consistent
            %so far. May be due to ROC calculation somewhere inside the code
            X =[X;ones(N-length(X),1)];
            Y = [Y; ones(N-length(Y),1)];
        else
            X = X(1:N);
            X(end) = 1.0;

            Y = Y(1:N);
            Y(end) = 1.0;
        end
        
        X_All{:,i} = X(:,1);
        Y_All{:,i} = Y(:,1);
        AUC_All(i) = AUC(1);
    end

    X_mean = mean(cell2mat(X_All),2); %X_std = std(New_X_All,0,2);
    Y_mean = mean(cell2mat(Y_All),2); Y_std = std(cell2mat(Y_All),0,2);
    AUC_mean = mean(AUC_All,2);
    std(AUC_All)
    plot(X_mean,Y_mean, "Color",paramColor,'LineWidth',2.0);
    hold on
    patch([X_mean; flipud(X_mean)], [Y_mean+Y_std; flipud(Y_mean-Y_std)], paramColor, 'EdgeColor','none', 'FaceAlpha',0.2, 'HandleVisibility','off')
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