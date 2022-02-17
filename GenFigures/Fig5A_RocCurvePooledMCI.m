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

%% Aggregating the MCI sample
AggregateMCI;

%% Calculating our model
%no consider of the rotation gain factor in leg2, i.e., gamma, G3=1, g2=1, g3, k3, sigma, nu. #params=6
config.UseGlobalSearch = true;
config.ModelName = "BaseModel";
config.NumParams = 5;
[MCIParameters, ~, ~, ~, AllYoungIC_Base] = getResultsAllConditions(YoungControls, config);

%% A function for getting Results from All Conditions
function [AllParams, AllX, AllDX, AllTheta, AllIC] = getResultsAllConditions(TransformedData, config)
    %get the estimated parameters, X, DX, Theta, IC for all trial
    %conditions for each group of data.
    AllParams = cell(0); AllX = cell(0);AllDX = cell(0); AllTheta = cell(0); AllIC = cell(0);
    for TRIAL_FILTER=1:3
        config.TrialFilter = TRIAL_FILTER;
        tic
        disp('%%%%%%%%%%%%%%% PERFORMING FITTING %%%%%%%%%%%%%%%');
        Results = PerformGroupFit(TransformedData, config);

        AllParams{TRIAL_FILTER} = Results.estimatedParams; 
        AllX{TRIAL_FILTER}      = Results.X;
        AllDX{TRIAL_FILTER}     = Results.DX;      
        AllTheta{TRIAL_FILTER}  = Results.THETADX;
        AllIC{TRIAL_FILTER}     = Results.IC;
        
        toc
    end
end    