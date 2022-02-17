%% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
HealthyControls = TransformPaths(HealthyControls);
MCIPos          = TransformPaths(MCIPos);
MCINeg          = TransformPaths(MCINeg);
Unknown         = TransformPaths(Unknown);
savefolder = pwd + "/Output/";

%% setting the configuration
resultfolder = savefolder+"PaperFigs/Fig5A";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

ColorPattern;

%% Aggregating the MCI sample
AggregateMCI;

%% Calculating our model
%no consider of the rotation gain factor in leg2, i.e., gamma, G3=1, g2=1, g3, k3, sigma, nu. #params=6
%"VariableNames", {'Gamma','g','b','sigma','nu'});
config.UseGlobalSearch = true;
config.ModelName = "BaseModel";
config.NumParams = 5;
[MCIAllParameters, ~, ~, ~, ~] = getResultsAllConditions(MCIAll, config);
MCIAllParameters = MCIAllParameters(:,[1 4 5 6 7]);

[HealthyControlsParameters, ~, ~, ~, ~] = getResultsAllConditions(HealthyControls, config);
HealthyControlsParameters = HealthyControlsParameters(:,[1 4 5 6 7]);

%% Preparing the logistic regression
AllData = [HealthyControlsParameters; MCIAllParameters];
logicalResponse = (1:height(HealthyControlsParameters) + height(MCIAllParameters))' > height(HealthyControlsParameters);
hcLabel     = cell(height(HealthyControlsParameters),1);
hcLabel(:)  = {'HC'};
mciLabel    = cell(height(MCIAllParameters),1);
mciLabel(:) = {'MCI'};
labels      = [hcLabel;mciLabel];

clear hcLabel mciLabel

%% Fitting the logistic regression
mdl = fitglm(AllData,logicalResponse,'Distribution', 'binomial','Link','logit');

scores = mdl.Fitted.Probability;
[X,Y,T,AUC] = perfcurve(labels,scores,'MCI');

plot(X,Y);

%% A function for getting Results from All Conditions
function [AllParams, AllX, AllDX, AllTheta, AllIC] = getResultsAllConditions(TransformedData, config)
    %get the estimated parameters, X, DX, Theta, IC for all trial
    %conditions for each group of data.    
    config.TrialFilter = 0;
    tic
    disp('%%%%%%%%%%%%%%% PERFORMING FITTING %%%%%%%%%%%%%%%');
    Results = PerformGroupFit(TransformedData, config);
    toc

    AllParams = Results.estimatedParams;
    AllX      = Results.X;
    AllDX     = Results.DX;
    AllTheta  = Results.THETADX;
    AllIC     = Results.IC;
    
end    