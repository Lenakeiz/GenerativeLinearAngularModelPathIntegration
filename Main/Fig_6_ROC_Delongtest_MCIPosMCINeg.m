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

%% Collecting behavioral information from output
% To allow fair comparisons wer are going to remove participants for which
% our model did not fit. See Methods for details on participants were excluded from the analysis.
MCINegDistanceError          = averageAcrossConditions_Data(MCINeg.Results.PropDistErr, MCINeg.Results.estimatedParams);
MCIPosDistanceError          = averageAcrossConditions_Data(MCIPos.Results.PropDistErr, MCIPos.Results.estimatedParams);

MCINegAngErr                 = averageAcrossConditions_Data(MCINeg.Results.PropAngErr, MCINeg.Results.estimatedParams);
MCIPosAngErr                 = averageAcrossConditions_Data(MCIPos.Results.PropAngErr, MCIPos.Results.estimatedParams);

allBehavMCIPos    = [MCIPosDistanceError MCIPosAngErr];
allBehavMCINeg    = [MCINegDistanceError MCINegAngErr];

%% Collecting model information from output
MCINegParameters          = averageAcrossConditions_Model(MCINeg.Results.estimatedParams);
MCIPosParameters          = averageAcrossConditions_Model(MCIPos.Results.estimatedParams);

%%
disp("%%%%%%%%%%%%%%% ROC MCI+ vs MCI- - model parameters estimation %%%%%%%%%%%%%%%")
rng("default");

param_idx = 7;
Model_dat = [MCIPosParameters(:,param_idx);MCINegParameters(:,param_idx)];
label1     = cell(height(MCIPosParameters(:,param_idx)),1);label1(:)  = {'MCIPos'};
label2     = cell(height(MCINegParameters(:,param_idx)),1);label2(:)  = {'MCINeg'};
Model_label = [label1;label2];

data_idx = 1; %1 for linear; 2 for angular
Behav_dat = [allBehavMCIPos(:,data_idx);allBehavMCINeg(:,data_idx)];
label1     = cell(height(allBehavMCIPos(:,data_idx)),1);label1(:)  = {'MCIPos'};
label2     = cell(height(allBehavMCINeg(:,data_idx)),1);label2(:)  = {'MCINeg'};
Behav_label  = [label1;label2];

Delong_output = DelongtestwithSVM(Model_dat, Model_label, Behav_dat, Behav_label, 'MCIPos', 'MCINeg');

%carry out a two sample Kolmogorov-Smirnov test
% 
% norm_vector1 = (Delong_output.AUCbehav - mean(Delong_output.AUCbehav)) / std(Delong_output.AUCbehav);
% norm_vector2 = (Delong_output.AUCmodel - mean(Delong_output.AUCmodel)) / std(Delong_output.AUCmodel);
% [h, p, ksstat] = kstest2(norm_vector1, norm_vector2);


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

    postp_model = zeros(23,1000);
    postp_behav = zeros(23,1000);

    warning('off');

    for i_cross_valid = 1:M
        %model
        mdl_model  = fitcsvm(Model_dat,Model_label,ClassNames=[{param1Label},{param2Label}], Standardize=true, CrossVal="on", Holdout=0.4, KernelFunction="linear");
 
        trainingIDs = training(mdl_model.Partition);
        Xtraining = Model_dat(trainingIDs,:);
        Ytraining = Model_label(trainingIDs);
        
        testingIDs = test(mdl_model.Partition);
        Xtest = Model_dat(testingIDs,:);
        Ytest_model = Model_label(testingIDs);
        compactModel = mdl_model.Trained{1};

        comp_mdl_post = fitPosterior(compactModel,Xtraining,Ytraining);
        [~,post_probabilities_model] = predict(comp_mdl_post,Xtest);
        
        [X,Y,~,AUC] = perfcurve(Ytest_model, post_probabilities_model(:,2), {param2Label});
        AUC_All_model(i_cross_valid) = AUC(1);
        N = length(Ytest_model);
        postp_model(testingIDs,i_cross_valid) = post_probabilities_model(:,1);
     

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
        Xtraining = Behav_dat(trainingIDs,:);
        Ytraining = Behav_label(trainingIDs);
        Xtest = Behav_dat(testingIDs,:);
        Ytest_behav = Behav_label(testingIDs);

        mdl_behav  = fitcsvm(Behav_dat,Behav_label,ClassNames=[{param1Label},{param2Label}],Standardize=true, KernelFunction="linear");

        comp_mdl_post = fitPosterior(mdl_behav,Xtraining,Ytraining);
        
        [~,post_probabilities_data] = predict(comp_mdl_post,Xtest);
        [X,Y,~,AUC] = perfcurve(Ytest_behav, post_probabilities_data(:,2), {param2Label});
        N = length(Ytest_behav);
        postp_behav(testingIDs,i_cross_valid) = post_probabilities_data(:,1);

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
        countMCIPos = sum(cellfun(@(x) strcmp(x, 'MCIPos'), Ytest_behav));
        countMCINeg = sum(cellfun(@(x) strcmp(x, 'MCINeg'), Ytest_behav));
        sample.spsizes = [countMCIPos, countMCINeg];
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
    nonnanP = allP(~isnan(allP));
    p005ratio = sum(nonnanP<0.05)/length(nonnanP);

    Delong_output.pvalueratio = p005ratio;
    Delong_output.AUCmodel = AUC_All_model;
    Delong_output.AUCbehav = AUC_All_behav;

    %AUC on mean posterior values
    mean_post_model = mean(postp_model,2);
    %[X,Y,~,AUC_model] = perfcurve(Model_label, mean_post_model, {param1Label});

    mean_post_behav = mean(postp_behav,2);
    %[X,Y,~,AUC_behav] = perfcurve(Model_label, mean_post_behav, {param1Label});

    countMCIPos = sum(cellfun(@(x) strcmp(x, 'MCIPos'), Model_label));
    countMCINeg = sum(cellfun(@(x) strcmp(x, 'MCINeg'), Model_label));
    sample.spsizes = [countMCIPos, countMCINeg];
    sample.ratings = [mean_post_model';mean_post_behav'];

    [aucs, delongcov] = fastDeLong(sample);
    meanpost_pvalue = calpvalue(aucs, delongcov);
    Delong_output.meanpost_pvalue = meanpost_pvalue;
    Delong_output.delongcov = delongcov;
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