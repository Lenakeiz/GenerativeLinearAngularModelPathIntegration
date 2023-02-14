%% Preparing the data
VAM_PrepareBaseConfig

% Preprocessing the data
VAM_PreprocessData

% Setting the model we are interested in
% Eventually modify config paramteters we are interested in. For example
% for this graph we are not interested in running the model with splitted
% conditions so we will set the relative config to false 
% force it to not run
rng("default");
config.useTrialFilter = true;
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName); % Set 100 here to avoid producing the model
% Run the model
VAM

%% Model run completed, preparing the data for plotting figures
config.ResultFolder = pwd + "/Output/ModelFigures/Fig5/ModelParameters/All_SVM_Models/Grouped_beta_k_g2_g3_sigma_nu";
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

% Loading previously saved SVM
if(exist('Data/SVM_all.mat','file') == 2)
    load Data/SVM_all.mat
end
%% Generating the SVM models. Previous run will be saved into Data folder in SVM_all.mat
NSamples = 100;
KernelFunction = "linear";
disp("%%%%%%%%%%%%%%% ALL SVM MODELS %%%%%%%%%%%%%%%")
disp("%%%%%%%%%%%%%%% ROC MCI pooled vs HC - model parameters estimation %%%%%%%%%%%%%%%")
rng("default");
HCvsMCISVMModels = getSVMBestModels(HealthyControlsParameters,MCIAllParameters,'HC', 'MCI', NSamples, KernelFunction);
%plotSVMResults(HealthyControlsParameters,MCIAllParameters,'HC', 'MCI',config, plotInfo);
%
plotInfo.Title = "MCI - / MCI +";
disp("%%%%%%%%%%%%%%% ROC MCI positive vs MCI negative - model parameters estimation %%%%%%%%%%%%%%%")
rng("default");
MCIPosvsMCINegSVMModels = getSVMBestModels(MCINegParameters, MCIPosParameters,'MCIneg', 'MCIpos', NSamples, KernelFunction);

save Data/SVM_all.mat HCvsMCISVMModels MCIPosvsMCINegSVMModels
%% Plotting the SVM models in instagram. Each one for the parameter of the best model.

% Plotting variables
close all;
plotInfo.figurePosition = [100 100 750 150;...
                           100 100 750 180;...
                           100 100 750 180;...
                           100 100 750 180;...
                           100 100 750 180;...
                           100 100 450 125;...
                           ];
plotInfo.defaultTextSize = 18;
plotInfo.defaultLineSize = 1.8;
plotInfo.titleFontSize = 16;
plotInfo.XlabelSize = [16 16 16 16 16 16];
plotInfo.YlabelSize = [14 14 14 14 14 14];
plotInfo.XaxisSize  = [13 13 13 13 12 13];
plotInfo.YaxisSize  = [13 13 13 13 13 13];
plotInfo.lineAlpha = 0.6;
plotInfo.textOffset = 0.1;
plotInfo.TextSize = [10 10 8 10 10 10];
plotInfo.yTickValues = [{0:0.5:1};...
                  {0:0.5:1};...
                  {0:0.5:1};...
                  {0:0.5:1};...
                  {0:0.5:1};...
                  {0:0.5:1};...
    ];
plotInfo.XlabelRotation = [90 90 90 90 90 0];
plotInfo.YLabel = "True positive rate";
plotInfo.XLabel = "False positive rate";
plotInfo.Title = "Healthy controls / pooled MCI";
plotInfo.visible = "on";

plotInfo.Title = "HC / MCI";
plotInfo.barColor = config.color_scheme_npg(5,:);
plotSVMResults(HealthyControlsParameters, MCIAllParameters,'HC', 'MCI', config, plotInfo, HCvsMCISVMModels);
%
plotInfo.Title = "MCIneg / MCIpos";
plotInfo.barColor = config.color_scheme_npg(6,:);
plotSVMResults(MCINegParameters, MCIPosParameters,'MCIneg', 'MCIpos', config, plotInfo, MCIPosvsMCINegSVMModels);
%%
function FS_out = getSVMBestModels(param1, param2, param1Label, param2Label, NSamples, kernelFunction)
    
    % Preparing data for the svm input
    allData = [param1; param2];
    allDatalogicalResponse = (1:height(param1) + height(param2))' > height(param1);
    label1     = cell(height(param1),1);
    label1(:)  = {param1Label};
    label2     = cell(height(param2),1);
    label2(:)  = {param2Label};
    allLabels  = [label1;label2];
    
    clear param1 param2

    % set the holdout on training/testing for cross validation
    holdout = 0.40;

    n_params = size(allData,2);
    params_idx = 1:n_params;

    FS_out = cell(n_params,1);

    % We look through combination of parameters starting from one
    for combination_idx = 1:n_params
        params_combination = nchoosek(params_idx,combination_idx);

        FS_out{combination_idx} = struct;
        FS_out{combination_idx}.n_params = combination_idx;
        FS_out{combination_idx}.n_combinations = height(params_combination);
        FS_out{combination_idx}.models = cell(height(params_combination),1);
        % from the current set of grouped parameters let s find the one
        % that gives the best AUC
        for curr_comb_idx = 1:height(params_combination)
                      
            FS_out{combination_idx}.models{curr_comb_idx} = struct;
            FS_out{combination_idx}.models{curr_comb_idx}.AUCmean = 0;
            FS_out{combination_idx}.models{curr_comb_idx}.AUCstd = 0;
            FS_out{combination_idx}.models{curr_comb_idx}.indeces = 0;
            FS_out{combination_idx}.models{curr_comb_idx}.X_mean = 0;
            FS_out{combination_idx}.models{curr_comb_idx}.Y_mean = 0;
            FS_out{combination_idx}.models{curr_comb_idx}.Y_std = 0;
          
            M = NSamples;
            X_All = cell(1,M);
            Y_All = cell(1,M);
            AUC_All = zeros(1,M);
            
            curr_filtered_all_data = allData(:,params_combination(curr_comb_idx,:));

            warning('off');
            for i_cross_valid = 1:M

                % fitting 
                mdl  = fitcsvm(curr_filtered_all_data,allLabels,ClassNames=[{param1Label},{param2Label}],Standardize=true, CrossVal="on", Holdout=holdout, KernelFunction=kernelFunction);
                trainingIDs = training(mdl.Partition);
                Xtraining = curr_filtered_all_data(trainingIDs,:);
                Ytraining = allLabels(trainingIDs);
                testingIDs = test(mdl.Partition);
                Xtest = curr_filtered_all_data(testingIDs,:);
                Ytest = allLabels(testingIDs);

                clear trainingIDs testingIDs
                compactModel = mdl.Trained{1};

                comp_mdl_post = fitPosterior(compactModel,Xtraining,Ytraining);
                %predict on leave-out testing data
                [predicted_labels,post_probabilities] = predict(comp_mdl_post,Xtest);
                [X,Y,~,AUC] = perfcurve(Ytest, post_probabilities(:,2), {param2Label});
                AUC_All(i_cross_valid) = AUC(1);
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

                X_All{:,i_cross_valid} = X(:,1);
                Y_All{:,i_cross_valid} = Y(:,1);

            end

            warning("on");
            X_mean = mean(cell2mat(X_All),2,"omitnan");
            Y_mean = mean(cell2mat(Y_All),2,"omitnan");
            Y_std = std(cell2mat(Y_All),0,2,"omitnan")./sqrt(N);

            FS_out{combination_idx}.models{curr_comb_idx}.AUCmean = mean(AUC_All,2,"omitnan");
            FS_out{combination_idx}.models{curr_comb_idx}.AUCstd = std(AUC_All,0,"omitnan");
            FS_out{combination_idx}.models{curr_comb_idx}.indeces = params_combination(curr_comb_idx,:);
            FS_out{combination_idx}.models{curr_comb_idx}.X_mean = X_mean;
            FS_out{combination_idx}.models{curr_comb_idx}.Y_mean = Y_mean;
            FS_out{combination_idx}.models{curr_comb_idx}.Y_std = Y_std;

            disp(['Fitting SVM on params ',...
                sprintf('%d,',params_combination(curr_comb_idx,:)),...
                ' size: ',...
                sprintf('%d,', size(curr_filtered_all_data)),...
                ' AUC = ',...
                sprintf('%.3f,', FS_out{combination_idx}.models{curr_comb_idx}.AUCmean),...
                ]);

        end
        clear curr_comb_idx

    end
    clear n_params

    warning('on')

end

%%
function plotSVMResults(params1, params2, params1groupName, params2groupName, config, plotInfo, SVMModels)
    %
    parametersName = [{'\beta'}, {'\itk'}, {'g_2'}, {'g_3'}, {'\sigma'}, {'\nu'}];
    %plotting a figure for each of the models
    for i_param =  1:length(config.ParamName)

        % set figure info
        %f = figure('visible','off','Position', [100 100 1000 500]);
        f = figure('visible', plotInfo.visible, 'Position', plotInfo.figurePosition(i_param,:));
        %%% Font type and size setting %%%
        % Using Arial as default because all journals normally require the font to
        % be either Arial or Helvetica
        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',plotInfo.axisSize)
        set(0,'DefaultTextFontSize',plotInfo.axisSize)
        
        hold on

        barValues = cellfun(@(x) x.AUCmean, SVMModels{i_param}.models);
        barSem = cellfun(@(x) x.AUCstd, SVMModels{i_param}.models);
        
        b = bar(barValues);
        eb = errorbar(1:numel(barValues), barValues, barSem, '.', 'color', 'k', 'LineWidth', 2);

        b.FaceColor = plotInfo.barColor;
        
        for i = 1:numel(barValues)
            xpos = b.XData(i); % + b.BarWidth/2;
            ypos = barValues(i) + barSem(i) + plotInfo.textOffset;
            text(xpos, ypos, num2str(round(barValues(i),2), '%.2f'), HorizontalAlignment = 'center', FontSize=plotInfo.TextSize(i_param));
        end

        hold off

        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.1 .1 .1], ...
            'YColor'      , [.1 .1 .1] );

        ax = gca;
        ax.XTick  = 1:SVMModels{i_param}.n_combinations;
        
        xlabelstringsindices = cellfun(@(x) x.indeces, SVMModels{i_param}.models, UniformOutput=false);
        labels = cell(1,height(xlabelstringsindices));
        for i_label = 1:height(xlabelstringsindices)
            labels{1,i_label} = getLabelFromIndeces(xlabelstringsindices{i_label}(1,:), parametersName);
        end
    
        ax.XTickLabel = labels;
        ax.XTickLabelRotation = plotInfo.XlabelRotation(i_param);

        ax.YTick = plotInfo.yTickValues{i_param}(:);
        ax.YLabel.String = 'AUC';
        ax.LineWidth = plotInfo.defaultLineSize;
        ax.XAxis.FontSize = plotInfo.XaxisSize(i_param);
        ax.YAxis.FontSize = plotInfo.YaxisSize(i_param);
        ax.XLabel.FontSize = plotInfo.XlabelSize(i_param);
        ax.YLabel.FontSize = plotInfo.YlabelSize(i_param);

        drawnow;

        exportName = [params1groupName 'vs' params2groupName '_' num2str(i_param) 'params'];

        exportgraphics(f,config.ResultFolder+"/"+convertCharsToStrings(exportName)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/"+convertCharsToStrings(exportName)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end

%% Extract the label from the indeces
function labelOut = getLabelFromIndeces(indecesArray, parametersName)
    labelOut = '';
    for i_text = 1:length(indecesArray)
        labelOut = [labelOut ' ' parametersName{indecesArray(i_text)}];
    end
    labelOut = strtrim(labelOut);
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