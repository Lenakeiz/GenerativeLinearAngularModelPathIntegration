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

%% Collecting behavioral information from output
% To allow fair comparisons wer are going to remove participants for which
% our model did not fit. See Methods for details on participants were excluded from the analysis.

HealthyControlsDistanceError = averageAcrossConditions_Data(HealthyControls.Results.PropDistErr, HealthyControls.Results.estimatedParams);
MCIUnkDistanceError          = averageAcrossConditions_Data(MCIUnk.Results.PropDistErr, MCIUnk.Results.estimatedParams);
MCINegDistanceError          = averageAcrossConditions_Data(MCINeg.Results.PropDistErr, MCINeg.Results.estimatedParams);
MCIPosDistanceError          = averageAcrossConditions_Data(MCIPos.Results.PropDistErr, MCIPos.Results.estimatedParams);

MCIAllDistanceError          = [MCIUnkDistanceError; MCINegDistanceError; MCIPosDistanceError];

HealthyControlsAngErr        = averageAcrossConditions_Data(HealthyControls.Results.PropAngErr, HealthyControls.Results.estimatedParams);
MCIUnkAngErr                 = averageAcrossConditions_Data(MCIUnk.Results.PropAngErr, MCIUnk.Results.estimatedParams);
MCINegAngErr                 = averageAcrossConditions_Data(MCINeg.Results.PropAngErr, MCINeg.Results.estimatedParams);
MCIPosAngErr                 = averageAcrossConditions_Data(MCIPos.Results.PropAngErr, MCIPos.Results.estimatedParams);
MCIAllAngErr                 = [MCIUnkAngErr; MCINegAngErr; MCIPosAngErr];

allBehavHC        = [HealthyControlsDistanceError HealthyControlsAngErr];
allBehavPooledMCI = [MCIAllDistanceError MCIAllAngErr];

%% Collecting model information from output
HealthyControlsParameters = averageAcrossConditions_Model(HealthyControls.Results.estimatedParams);
MCIUnkParameters          = averageAcrossConditions_Model(MCIUnk.Results.estimatedParams);
MCINegParameters          = averageAcrossConditions_Model(MCINeg.Results.estimatedParams);
MCIPosParameters          = averageAcrossConditions_Model(MCIPos.Results.estimatedParams);
MCIAllParameters          = [MCIUnkParameters; MCINegParameters; MCIPosParameters]; 

%%
disp("%%%%%%%%%%%%%%% ROC MCI pooled vs HC - model parameters estimation %%%%%%%%%%%%%%%")
rng("default");

param_idx = 1;
Model_dat = [HealthyControlsParameters(:,param_idx);MCIAllParameters(:,param_idx)];
label1     = cell(height(HealthyControlsParameters(:,param_idx)),1);label1(:)  = {'HC'};
label2     = cell(height(MCIAllParameters(:,param_idx)),1);label2(:)  = {'MCI'};
Model_label = [label1;label2];

data_idx = 2; %1 for linear; 2 for angular
Behav_dat = [allBehavHC(:,data_idx);allBehavPooledMCI(:,data_idx)];
label1     = cell(height(allBehavHC(:,data_idx)),1);label1(:)  = {'HC'};
label2     = cell(height(allBehavPooledMCI(:,data_idx)),1);label2(:)  = {'MCI'};
Behav_label  = [label1;label2];

HCvsMCISVMModels_SingleParameters = DelongtestwithSVM(Model_dat, Model_label, Behav_dat, Behav_label, 'HC', 'MCI');


%%
function Delong_output = DelongtestwithSVM(Model_dat, Model_label, Behav_dat, Behav_label, param1Label, param2Label)
    
    M = 1000;
    X_All_model = cell(1,M);
    Y_All_model = cell(1,M);
    AUC_All_model = zeros(1,M);

    X_All_behav = cell(1,M);
    Y_All_behav = cell(1,M);
    AUC_All_behav = zeros(1,M);

    allPvalue = cell(1,M);

    warning('off');

    for i_cross_valid = 1:M
        %model
        mdl_model  = fitcsvm(Model_dat,Model_label,ClassNames=[{param1Label},{param2Label}], Standardize=true, CrossVal="on", Holdout=0.4, KernelFunction="linear");
 
        trainingIDs = training(mdl_model.Partition);
        Xtraining = Model_dat(trainingIDs,:);
        Ytraining = Model_label(trainingIDs);
        
        testingIDs = test(mdl_model.Partition);
        Xtest = Model_dat(testingIDs,:);
        Ytest = Model_label(testingIDs);
        compactModel = mdl_model.Trained{1};

        comp_mdl_post = fitPosterior(compactModel,Xtraining,Ytraining);
        [~,post_probabilities_model] = predict(comp_mdl_post,Xtest);
        
        [X,Y,~,AUC] = perfcurve(Ytest, post_probabilities_model(:,2), {param2Label});
        AUC_All_model(i_cross_valid) = AUC(1);
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

        X_All_model{:,i_cross_valid} = X(:,1);
        Y_All_model{:,i_cross_valid} = Y(:,1);
        
        %data
        mdl_behav  = fitcsvm(Behav_dat,Behav_label,ClassNames=[{param1Label},{param2Label}],Standardize=true, CrossVal="on", Holdout=0.4, KernelFunction="linear");
        trainingIDs = training(mdl_behav.Partition);
        Xtraining = Behav_dat(trainingIDs,:);
        Ytraining = Behav_label(trainingIDs);

        testingIDs = test(mdl_behav.Partition);
        Xtest = Behav_dat(testingIDs,:);
        Ytest = Behav_label(testingIDs);
        compactModel = mdl_behav.Trained{1};

        comp_mdl_post = fitPosterior(compactModel,Xtraining,Ytraining);
        
        [~,post_probabilities_data] = predict(comp_mdl_post,Xtest);
        [X,Y,~,AUC] = perfcurve(Ytest, post_probabilities_data(:,2), {param2Label});
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

        X_All_behav{:,i_cross_valid} = X(:,1);
        Y_All_behav{:,i_cross_valid} = Y(:,1);

        AUC_All_behav(i_cross_valid) = AUC(1);
        
        %Delong test
        %make the sample
        countHC = sum(cellfun(@(x) strcmp(x, 'HC'), Ytest));
        countMCI = sum(cellfun(@(x) strcmp(x, 'MCI'), Ytest));
        sample.spsizes = [countHC, countMCI];
        sample.ratings = [post_probabilities_model(:,1)';post_probabilities_data(:,1)'];

        [aucs, delongcov] = fastDeLong(sample);
        pvalue = calpvalue(aucs, delongcov);
        allPvalue{:, i_cross_valid} = pvalue;

    end

    warning("on");

    X_mean_model = mean(cell2mat(X_All_model),2,"omitnan");
    Y_mean_model = mean(cell2mat(Y_All_model),2,"omitnan");
    Y_std_model = std(cell2mat(Y_All_model),0,2,"omitnan")./sqrt(N);
    AUCOut_model.mean = mean(AUC_All_model,2,"omitnan");
    AUCOut_model.CI   = [0 0];
    AUCOut_model.std  = std(AUC_All_model,0,"omitnan");

    X_mean_behav = mean(cell2mat(X_All_behav),2,"omitnan");
    Y_mean_behav = mean(cell2mat(Y_All_behav),2,"omitnan");
    Y_std_behav = std(cell2mat(Y_All_behav),0,2,"omitnan")./sqrt(N);
    AUCOut_behav.mean = mean(AUC_All_behav,2,"omitnan");
    AUCOut_behav.CI   = [0 0];
    AUCOut_behav.std  = std(AUC_All_behav,0,"omitnan");

    allP = cell2mat(allPvalue);
    p005 = sum(allP<0.05);

end





%% 
function dataout = averageAcrossConditions_Model(data)
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

%%
function dataout = averageAcrossConditions_Data(data, glampi_data)
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