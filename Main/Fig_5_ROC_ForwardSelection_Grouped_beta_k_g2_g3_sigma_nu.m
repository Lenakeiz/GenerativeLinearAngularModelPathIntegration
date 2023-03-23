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
config.ResultFolder = pwd + "/Output/ModelFigures/Fig5/ModelParameters/ForwardSelection/Grouped_beta_k_g2_g3_sigma_nu";
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
if(exist('Data/SVM_ForwardSelection.mat','file') == 2)
    load Data/SVM_ForwardSelection.mat
end

%% Plotting roc curve HC vs pooled MCI and MCI negative vs MCI positive

NSamples= 1000;
disp("%%%%%%%%%%%%%%% SVM FORWARD SELECTION %%%%%%%%%%%%%%%")
disp("%%%%%%%%%%%%%%% ROC MCI pooled vs HC - model parameters estimation %%%%%%%%%%%%%%%")
rng("default");
HCvsMCISVMModels_FS = getSVMBestModels(HealthyControlsParameters,MCIAllParameters,'HC', 'MCI', NSamples);
disp("%%%%%%%%%%%%%%% ROC MCI positive vs MCI negative - model parameters estimation %%%%%%%%%%%%%%%")
rng("default");
MCIPosvsMCINegSVMModels_FS = getSVMBestModels(MCINegParameters, MCIPosParameters,'MCIneg', 'MCIpos', NSamples);
clear NSamples
save Data/SVM_ForwardSelection.mat HCvsMCISVMModels_FS MCIPosvsMCINegSVMModels_FS
%%
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
plotInfo.Title = "Healthy controls / pooled MCI";
plotInfo.visible = "on";
plotInfo.figurePosition = [200 200 250 200];
%
plotInfo.Title = "HC / MCI";
plotSVMResults(HealthyControlsParameters,MCIAllParameters,'HC', 'MCI',config, plotInfo, HCvsMCISVMModels_FS);
%
plotInfo.Title = "MCI - / MCI +";
plotSVMResults(MCINegParameters, MCIPosParameters,'MCIneg', 'MCIpos', config, plotInfo, MCIPosvsMCINegSVMModels_FS);  
%%
function FS_out = getSVMBestModels(param1, param2, param1Label, param2Label, NSamples)
    
    % pre allocating output variable
    FS_out=struct;
    % taking track of the 4 best parameters
    FS_out.models = cell(4,1);
    for i = 1:4
        FS_out.models{i,1}.mean = 0;
        FS_out.models{i,1}.CI = [0 0];
        FS_out.models{i,1}.std = 0;
        % tracking the best model out of the forward search
        FS_out.models{i,1}.indeces = nan;
        FS_out.models{i,1}.X_mean = nan;
        FS_out.models{i,1}.Y_mean = nan;
        FS_out.models{i,1}.Y_std = nan;
    end
    
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

    fwd_search = struct();
    fwd_search.selected_from_previous = nan;
    
    % We look through combination of parameters starting from one
    for combination_idx = 1:n_params
        params_combination = nchoosek(params_idx,combination_idx);

        %From previous iteration filtering the rows that have the parameters with
        %best AUC
        if(~isnan(fwd_search.selected_from_previous))
            row_with_previous_found_indeces = any(ismember(params_combination, fwd_search.selected_from_previous), 2);
            params_combination = params_combination(row_with_previous_found_indeces,:);
        end
        clear row_with_previous_found_indeces

        % from the current set of grouped parameters let s find the one
        % that gives the best AUC
        for curr_comb_idx = 1:height(params_combination)
                
            M = NSamples;
            X_All = cell(1,M);
            Y_All = cell(1,M);
            AUC_All = zeros(1,M);
            
            curr_filtered_all_data = allData(:,params_combination(curr_comb_idx,:));

            warning('off');
            for i_cross_valid = 1:M

                % fitting 
                mdl  = fitcsvm(curr_filtered_all_data,allLabels,ClassNames=[{param1Label},{param2Label}],Standardize=true, CrossVal="on", Holdout=holdout, KernelFunction="linear");
                trainingIDs = training(mdl.Partition);
                Xtraining = curr_filtered_all_data(trainingIDs,:);
                Ytraining = allLabels(trainingIDs);
                testingIDs = test(mdl.Partition);
                Xtest = curr_filtered_all_data(testingIDs,:);
                Ytest = allLabels(testingIDs);

                clear trainingIDs testingIDs
                compactModel = mdl.Trained{1};

                comp_mdl_post = fitPosterior(compactModel,Xtraining,Ytraining);
                %comp_mdl_post = fitPosterior(compactModel,curr_filtered_all_data,allLabels);
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

            candidate_model{1,1}.mean = mean(AUC_All,2,"omitnan");
            candidate_model{1,1}.CI = [0 0];
            candidate_model{1,1}.std = std(AUC_All,0,"omitnan");
            candidate_model{1,1}.indeces = params_combination(curr_comb_idx,:);
            candidate_model{1,1}.X_mean = X_mean;
            candidate_model{1,1}.Y_mean = Y_mean;
            candidate_model{1,1}.Y_std = Y_std;

            disp(['Fitted SVM on params ',...
                sprintf('%d,',params_combination(curr_comb_idx,:)),...
                ' size: ',...
                sprintf('%d,', size(curr_filtered_all_data)),...
                ' AUC = ',...
                sprintf('%.3f,', candidate_model{1,1}.mean),...
            ]);

            % Updating the structure containing four best models
            for i = 1:size(FS_out.models,1)
                if(candidate_model{1,1}.mean > FS_out.models{i,1}.mean)
                    if(i == 1)
                        FS_out.models = [candidate_model;FS_out.models(i:end-1)];
                        % Only the best one will dominate the forward
                        % search model
                        fwd_search.selected_from_previous = params_combination(curr_comb_idx,:);
                    elseif(i == size(FS_out.models,1))
                        FS_out.models = [FS_out.models(1:end-1); candidate_model];
                    else
                        FS_out.models = [FS_out.models(1:i-1);candidate_model;FS_out.models(i:end-1)];
                    end
                    break;
                end
            end

        end
        clear curr_comb_idx

    end
    clear n_params

    warning('on')

end

%%
function plotSVMResults(params1, params2, params1groupName, params2groupName, config, plotInfo, SVMModels)
    %
    parametersName = [{'\beta'}, {'k'}, {'g_2'}, {'g_3'}, {'\sigma'}, {'\nu'}];
    colors = config.color_scheme_npg([8 3 6 9 2 4],:);
    
    % set figure info
    f = figure('visible', plotInfo.visible, 'Position', plotInfo.figurePosition);
    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
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

        model_parameters_name = parametersName(:,SVMModels.models{i,1}.indeces);
        l_text = "AUC(";
        for i_text = 1:size(model_parameters_name,2)
            l_text = l_text + convertCharsToStrings(model_parameters_name{i_text});
            if i_text<size(model_parameters_name,2)
                l_text = l_text + " ";
            end
        end
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
    %t= title(plotInfo.Title);
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
    
    clear parametersName filesName i colors f ll legendText found_parameters_name
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