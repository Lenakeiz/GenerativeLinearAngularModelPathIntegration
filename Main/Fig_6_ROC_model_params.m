%% Script to create output for Fig. 5 - receiver operating characteristic curve for each GLAMPI parameter 
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
rng("default");
config.useTrialFilter = true;
config.ModelName        =   "beta_k_g2_g3_m3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "m3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

% Preparing output
config.ResultFolder = pwd + "/Output/Fig6/ModelParameters/beta_k_g2_g3_m3_sigma_nu";
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Generating color scheme for our paper
ColorPattern;

%%
% Collecting information from output
YoungControlsParameters   = averageAcrossConditions(YoungControls.Results.estimatedParams);
HealthyControlsParameters = averageAcrossConditions(HealthyControls.Results.estimatedParams);
MCIUnkParameters          = averageAcrossConditions(MCIUnk.Results.estimatedParams);
MCINegParameters          = averageAcrossConditions(MCINeg.Results.estimatedParams);
MCIPosParameters          = averageAcrossConditions(MCIPos.Results.estimatedParams);

MCIAllParameters          = [MCIUnkParameters; MCINegParameters; MCIPosParameters]; 

% Loading previously saved SVM
if(exist('Data/ROC_SVM_ModelParameters_results.mat','file') == 2)
    load Data/ROC_SVM_ModelParameters_results.mat
end

% Samples used for cross validation
NSamples= 1000;
% Kernel for the SVM
KernelFunction = "linear";
Holdout = 0.4;
disp("%%%%%%%%%%%%%%% Fitting SVM on parameters %%%%%%%%%%%%%%%");
disp("%%%%%%%%%%%%%%% ROC MCI pooled vs HC - model parameters estimation %%%%%%%%%%%%%%%")
rng("default");
HCvsMCISVMModels_SingleParameters = getSVMModels(HealthyControlsParameters,MCIAllParameters,'HC', 'MCI', config,  NSamples, KernelFunction, Holdout, "HOMCI");
disp("%%%%%%%%%%%%%%% ROC MCI positive vs MCI negative - model parameters estimation %%%%%%%%%%%%%%%")
rng("default");
MCIPosvsMCINegSVMModels_SingleParameters = getSVMModels(MCINegParameters, MCIPosParameters,'MCIneg', 'MCIpos', config, NSamples, KernelFunction, Holdout, "MCI");
clear NSamples
save Data/ROC_SVM_ModelParameters_results.mat HCvsMCISVMModels_SingleParameters MCIPosvsMCINegSVMModels_SingleParameters
disp("%%%%%%%%%%%%%%% SVM fitting done %%%%%%%%%%%%%%%");

% Plot results
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
plotInfo.Title = "Healthy controls / pooled MCI";
plotInfo.visible = "off";
plotInfo.figurePosition = [200 200 250 200];

disp("%%%%%%%%%%%%%%% Plotting ROC HC vs MCI model %%%%%%%%%%%%%%%");
plotInfo.Title = "HC / MCI";
plotSVMResults(HealthyControlsParameters,MCIAllParameters,'HC', 'MCI', config, plotInfo, HCvsMCISVMModels_SingleParameters);

disp("%%%%%%%%%%%%%%% Plotting ROC MCI positive vs MCI negative %%%%%%%%%%%%%%%");
plotInfo.Title = "MCI - / MCI +";
plotSVMResults(MCINegParameters, MCIPosParameters,'MCIneg', 'MCIpos', config, plotInfo, MCIPosvsMCINegSVMModels_SingleParameters);  

% Final cleanup to leave workspace as the end of the Preprocessing stage.
% Remove if you want to take a look at the output data.
clearvars -except config YoungControls HealthyControls MCINeg MCIPos MCIUnk

%%
function SVM_out = getSVMModels(param1, param2, param1Label, param2Label, config, NSamples, kernelFunction, Holdout, type)
    
    % Preparing data for the svm input
    allData = [param1; param2];
    allDatalogicalResponse = (1:height(param1) + height(param2))' > height(param1);
    label1     = cell(height(param1),1);
    label1(:)  = {param1Label};
    label2     = cell(height(param2),1);
    label2(:)  = {param2Label};
    allLabels  = [label1;label2];
    holdout = Holdout;

    clear param1 param2

    n_params = size(allData,2);
    % Pre allocating output variable with empty values
    SVM_out=struct;
    SVM_out.models = cell(config.NumParams ,1);
    for i = 1:config.NumParams 
        SVM_out.models{i,1}.mean = 0;
        SVM_out.models{i,1}.std = 0;
        SVM_out.models{i,1}.index = nan;
        SVM_out.models{i,1}.X_mean = nan;
        SVM_out.models{i,1}.Y_mean = nan;
        SVM_out.models{i,1}.Y_std = nan;
    end

    % Looping through each of the parameters
    Delong_all = cell(7,1);
    for i_params = 1:n_params

        M = NSamples;
        X_All = cell(1,M);
        Y_All = cell(1,M);
        AUC_All = zeros(1,M);

        curr_filtered_all_data = allData(:,i_params);

        warning('off');

        % Cross-validation
        if type=="HOMCI"
            Ratings = zeros(1000,28);
        else
            Ratings = zeros(1000,9);
        end

        for i_cross_valid = 1:M

            mdl  = fitcsvm(curr_filtered_all_data,allLabels,ClassNames=[{param1Label},{param2Label}], Standardize=true, CrossVal="on", Holdout=holdout, KernelFunction=kernelFunction);
            trainingIDs = training(mdl.Partition);
            Xtraining = curr_filtered_all_data(trainingIDs,:);
            Ytraining = allLabels(trainingIDs);
            testingIDs = test(mdl.Partition);
            Xtest = curr_filtered_all_data(testingIDs,:);
            Ytest = allLabels(testingIDs);

            clear trainingIDs testingIDs
            compactModel = mdl.Trained{1};

            comp_mdl_post = fitPosterior(compactModel,Xtraining,Ytraining);
            [predicted_labels,post_probabilities] = predict(comp_mdl_post,Xtest);
            
            [X,Y,~,AUC, ~, subY] = perfcurve(Ytest, post_probabilities(:,2), {param2Label});
            AUC_All(i_cross_valid) = AUC(1);
            N = length(Ytest);

            Ratings(i_cross_valid,:) = post_probabilities(:,1)';

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

            X_All{:,i_cross_valid} = X(:,1);
            Y_All{:,i_cross_valid} = Y(:,1);

        end
        
        if type == "HOMCI"
            class1 = Ratings(:,1:13);
            class2 = Ratings(:,14:28);
        else
            class1 = Ratings(:,1:6);
            class2 = Ratings(:,7:9);
        end
        Delong_all{i_params} = [class1(:)', class2(:)'];

        warning("on");
        X_mean = mean(cell2mat(X_All),2,"omitnan");
        Y_mean = mean(cell2mat(Y_All),2,"omitnan");
        Y_std = std(cell2mat(Y_All),0,2,"omitnan")./sqrt(N);

        SVM_out.models{i_params,1}.mean = mean(AUC_All,2,"omitnan");
        SVM_out.models{i_params,1}.std = std(AUC_All,0,"omitnan");
        SVM_out.models{i_params,1}.index = config.ParamName{i_params};
        SVM_out.models{i_params,1}.X_mean = X_mean;
        SVM_out.models{i_params,1}.Y_mean = Y_mean;
        SVM_out.models{i_params,1}.Y_std = Y_std;

        disp(['Fitted SVM on params ',...
            sprintf('%s,',config.ParamName{i_params}),...
            ' AUC = ',...
            sprintf('%.3f,', SVM_out.models{i_params,1}.mean),...
            ]);
    end
    warning('on')
    
    Delong_all = cell2mat(Delong_all);
    %save for later Delong test
    if type == "HOMCI"
        sample.ratings = Delong_all;
        sample.spsizes = [13000,15000];
    else
        sample.ratings = Delong_all;
        sample.spsizes = [6000,3000]; 
    end
    save("Data/Delong_Model_"+type+".mat", "sample")


end

%%
function plotSVMResults(params1, params2, params1groupName, params2groupName, config, plotInfo, SVMModels)
    %
    parametersName = [{'\beta'}, {'k'}, {'g_2'}, {'g_3'}, {'m_3'}, {'\sigma'}, {'\nu'}];
    colors = config.color_scheme_npg([4 8 2 3 5 9 7],:);
    
    f = figure('visible', plotInfo.visible, 'Position', plotInfo.figurePosition);

    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',plotInfo.axisSize)
    set(0,'DefaultTextFontSize',plotInfo.axisSize)
    
    hold on

    for i =1:size(SVMModels.models,1)
        
        plot(SVMModels.models{i,1}.X_mean,SVMModels.models{i,1}.Y_mean,"Color",colors(i,:),'LineWidth',plotInfo.defaultLineSize);
        patch([SVMModels.models{i,1}.X_mean; flipud(SVMModels.models{i,1}.X_mean)],...
            [SVMModels.models{i,1}.Y_mean+SVMModels.models{i,1}.Y_std; flipud(SVMModels.models{i,1}.Y_mean-SVMModels.models{i,1}.Y_std)],...
            colors(i,:), 'EdgeColor','none', 'FaceAlpha',0.2, 'HandleVisibility','off');

        l_text = "AUC(";
        l_text = l_text + convertCharsToStrings(parametersName(i));
        clear i_text
        l_text = l_text + ") = " + num2str(round(SVMModels.models{i,1}.mean,2),'%.2f');
        legendText{1,i} = l_text;
        disp("AUC std : " + num2str(SVMModels.models{i,1}.std));
        clear l_text
    end
    clear i;
    
    hline = refline([1 0]);
    hline.Color = 'k';
    hline.LineWidth = 2;
    hline.LineStyle = '--';
    
    hold off
    
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
    
    clear parametersName filesName i colors f ll legendText found_parameters_name
end

%% 
function dataout = averageAcrossConditions(data)
    % average across the conditions, assuming a specific structure to this
    % data
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

    % remove nans if the row after merging data
    dataout = removeNanRows(dataout);

end