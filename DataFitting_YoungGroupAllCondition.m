%% FittingData 
% This script will lo%% FittingData 
% This script will load the data, rotate the paths and call the solver to find the parameters 
% that minimize the mean distance to generated x3' which matchs the physical leg 3

%% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
%load('Data/AllDataErrors2018_V2.mat');
%V3 including the missing participants, which should always being loaded now
load('Data/AllDataErrors2018_V3.mat');
savefolder = pwd + "/Output/";

%% choose the model 

%the full model with G1, G2, G3, g2, g3 k3, sigma, nu. #params=8
%Model_Name = "base_add_G3"; numParams=8;  

%the base model with G1, G2, G3=1, g2, g3, k3, sigma, nu. #params=7
%Model_Name = "base"; numParams=7;   

%consider k3=0 in the excution turning error of leg 3, i.e., G1, G2, G3=1, g2, g3, k3=0, sigma, nu. #params=6
%Model_Name = "set_k3_0"; numParams=6;   

%no consider of execution error i.e., G1, G2, G3=1, g2, g3=1, k3=0, sigma, nu. #params=5
%Model_Name = "set_g3_1_k3_0"; numParams=5; 

%no consider of the rotation gain factor in leg2, i.e., G1, G2, G3=1, g2=1, g3, k3, sigma, nu. #params=6
Model_Name = "set_g2_1"; numParams=6;  

%no consider of the length gain factor in leg2, i.e., G1, G2=1, G3=1, g2, g3, k3, sigma, nu. #params=6
%Model_Name = "set_bG2_1"; numParams=6;   

%no consider of the length gain factor in leg1, i.e., G1=1, G2, G3=1, g2, g3, k3, sigma, nu. #params=6
%Model_Name = "set_bG1_1"; numParams=6;   

%same length gain factor in leg1 & 2, i.e., G1=G2, g2, g3, k3, sigma, nu. #params=6
%Model_Name = "set_same_bG"; numParams=6; 

%same rotation gain factor in leg2 & 3, i.e., G1, G2, g2=g3, k3=0, sigma, nu.#params=5
%Model_Name = "set_same_sg"; numParams=5; 

%same length gain, same rotation gain, i.e., G1=G2, g2=g3, k3=0, sigma, nu.#params=4
%Model_Name = "set_same_Gg"; numParams=4; 

%% initial a empty cell for storing all configuration (for passing varibles to functions easliy) 
config.ModelName = Model_Name;
config.NumParams = numParams;

config.USEOoBTrials = true;
if config.USEOoBTrials==true
    resultfolder = savefolder+"YoungGroup/"+Model_Name+"/OoB/";
else
    resultfolder = savefolder+"YoungGroup/"+Model_Name+"/NoOoB/";
end
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end
config.ResultFolder = resultfolder;
%

%for storing all IC information
%create storing folder for trajectory if not exist
icfolder = resultfolder+"IC/";
if ~exist(icfolder, 'dir')
   mkdir(icfolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over three different conditions for Young participants
%load the data
if config.USEOoBTrials==true
    YoungControls = TransformPathsOoB(YoungControls);
else
    YoungControls = TransformPaths(YoungControls);
end
%get fitting results for all conditions
[AllYoungParams, AllYoungX, AllYoungDX, AllYoungTheta, AllYoungIC] = getResultsAllConditions(YoungControls, config);
%storing IC
save(icfolder+"ICYoung.mat", 'AllYoungIC');

%% a function for 
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
        AllX{TRIAL_FILTER} = Results.X;
        AllDX{TRIAL_FILTER}=Results.DX;      
        AllTheta{TRIAL_FILTER}=Results.THETADX;
        AllIC{TRIAL_FILTER}=Results.IC;
        toc
    end
end    
