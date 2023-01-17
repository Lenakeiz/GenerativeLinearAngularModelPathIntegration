%% Preparing the data
VAM_PrepareBaseConfig

%% Preprocessing the data
VAM_PreprocessData

%% Setting the model we are interested in
% Eventually modify config paramteters we are interested in. For example
% for this graph we are not interested in running the model with splitted
% conditions so we will set the relative config to false 
% force it to not run
rng("default")
config.useTrialFilter = true;
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);
% Run the model
VAM

%% Model run completed, preparing the data for plotting figures
config.ResultFolder = pwd + "/Output/MRILinearMixedEffectModel/StepwiseGLM/Beta_k_g2_g3_sigma_nu";
% Create storing folder for trajectory if not exist
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

ColorPattern;

%% Getting Information from results:
YoungControlsParameters   = averageAcrossConditions(YoungControls.Results.estimatedParams);
HealthyControlsParameters = averageAcrossConditions(HealthyControls.Results.estimatedParams);
MCIUnkParameters          = averageAcrossConditions(MCIUnk.Results.estimatedParams);
MCINegParameters          = averageAcrossConditions(MCINeg.Results.estimatedParams);
MCIPosParameters          = averageAcrossConditions(MCIPos.Results.estimatedParams);

%% Preparing linear mixed effect table
HealthyControls_lmetable = createLMETable(HealthyControls.MRI,HealthyControlsParameters);
MCIUnk_lmetable = createLMETable(MCIUnk.MRI,MCIUnkParameters);
MCINeg_lmetable = createLMETable(MCINeg.MRI,MCINegParameters);
MCIPos_lmetable = createLMETable(MCIPos.MRI,MCIPosParameters);

MRIModelParamsDataTable = [HealthyControls_lmetable; MCIUnk_lmetable; MCINeg_lmetable; MCIPos_lmetable];

%Removing nans (either from Nan parameters model or missing mri data)
filter = isnan(MRIModelParamsDataTable.beta) | isnan(MRIModelParamsDataTable.MCI);
MRIModelParamsDataTable = MRIModelParamsDataTable(~filter,:);
clear filter
%% 
% For each of the parameter we are trying to fit we calculate a stepwise
% general linear model to check what could be the best volumetric areas to
% predict that fitted value
clc;
% Using the Destrieux atlas from free surfer (https://www.sciencedirect.com/science/article/pii/S1053811910008542?via%3Dihub)
destVolumesIdx = [13 14 15 18 21:22 29 47:79];
removeThreshold = 0.05;

parameterName = "beta";
modelBeta = performStepwiseGLM(MRIModelParamsDataTable,parameterName,destVolumesIdx,removeThreshold);
%
parameterName = "k";
modelG2 = performStepwiseGLM(MRIModelParamsDataTable,parameterName,destVolumesIdx,removeThreshold);
%
parameterName = "g2";
modelG2 = performStepwiseGLM(MRIModelParamsDataTable,parameterName,destVolumesIdx,removeThreshold);
%
parameterName = "g3";
modelG3 = performStepwiseGLM(MRIModelParamsDataTable,parameterName,destVolumesIdx,removeThreshold);
%
parameterName = "sigma";
modelSigma = performStepwiseGLM(MRIModelParamsDataTable,parameterName,destVolumesIdx,removeThreshold);
%
parameterName = "nu";
modelNu = performStepwiseGLM(MRIModelParamsDataTable,parameterName,destVolumesIdx,removeThreshold);
%%
function modelOut = performStepwiseGLM(MRIModelParamsDataTable,parameterName,destVolumesIdx,removeThreshold)
    %mdl = stepwiseglm(MRIModelParamsDataTable,"linear","ResponseVar", parameterName, "PredictorVars ", destVolumesIdx,"Criterion","PRemove",0.1)
    modelOut = stepwiselm(MRIModelParamsDataTable,"linear","Upper","linear","PredictorVars",destVolumesIdx,"ResponseVar",parameterName,"Criterion","aic","PRemove",removeThreshold)
end

%% Creates a unique table between mri and model parameters
function dataout = createLMETable(MRI, GroupParams)
    GroupParamsTable = array2table(GroupParams, "VariableNames", {'beta', 'k', 'g2', 'g3', 'sigma', 'nu'});
    dataout = [GroupParamsTable,MRI];
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

end