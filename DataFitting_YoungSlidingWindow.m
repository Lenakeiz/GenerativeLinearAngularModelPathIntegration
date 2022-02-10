%% FittingData 
% This script will load the data, rotate the paths and call the solver to find the parameters 
% that minimize the mean distance to generated x3' which matchs the physical leg 3

%% Cleaning variables
clearvars; clear all; close all; clc;

%% initial a empty cell for storing all configuration
config = cell(0);  

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');

%% Decide whether to add back those Out-of-Boundary trials
% Remember to change the PLOT_FEEDBACK to 0 inside the function inside otherwise it will pause and
% display the paths
disp('%%%%%%%%%%%%%%% TRANSFORME DATA ... %%%%%%%%%%%%%%%');
resultfolder = "C:/Users/Zilong/Desktop/path integration model/Andrea's matlab code/PathIntegrationDataEstimationModelling/Output/OoB/";
%Note that I put "CalculateOoBErrors" into TransformPathsOoB to remove the repeat calling of "CalculateOoBErrors" here.
YoungControls        = TransformPathsOoB(YoungControls);

config.ResultFolder = resultfolder;
%% choose the model that we want to optimize

%the full model with G1, G2, G3, g2, g3 k3, sigma, nu. #params=8
%Model_Name = "base_add_G3"; numParams=8;  

%the base model with G1, G2, G3=1, g2, g3, k3, sigma, nu. #params=7
Model_Name = "base"; numParams=7;   

%consider k3=0 in the excution turning error of leg 3, i.e., G1, G2, G3=1, g2, g3, k3=0, sigma, nu. #params=6
%Model_Name = "set_k3_0"; numParams=6;   

%no consider of execution error i.e., G1, G2, G3=1, g2, g3=1, k3=0, sigma, nu. #params=5
%Model_Name = "set_g3_1_k3_0"; numParams=5; 

%no consider of the rotation gain factor in leg2, i.e., G1, G2, G3=1, g2=1, g3, k3, sigma, nu. #params=6
%Model_Name = "set_g2_1"; numParams=6;  

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

config.ModelName = Model_Name;
config.NumParams = numParams;

%% loop over three different conditions

%initial empty cells to store results
AllYoungParams= cell(0);

% 0 all trials, 1 no change, 2 no distal cues, 3 no optic flow
for TRIAL_FILTER=1:3
    config.TrialFilter = TRIAL_FILTER;
    % Performing data fitting
    % For reproducibility purposes setting the same random seed every time
    % rng('default');
    tic
    disp('%%%%%%%%%%%%%%% PERFORMING FIT YOUNG CONTROLS %%%%%%%%%%%%%%%');
    YoungParamValues= PerformGroupFitYoungSlidingWindow(YoungControls, 0.7, config);
    AllYoungParams{TRIAL_FILTER} = YoungParamValues;
    toc
end

%% visualizing how the value of estimated parameters changes along with the varying of leg length
PlotShadedSlidingWindowOfYoung(AllYoungParams, config);