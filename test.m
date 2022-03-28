%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default');

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
savefolder = pwd + "/Output/";

%% set up the configuration
config.Speed.alpha = 0.9;                               %Paramanter for running speed calculation
config.Speed.timeOffsetAfterFlagReach = 2;              %Time to track after flag reached in seconds 
config.Speed.smoothWindow = 10;                         % tracking rate should be 10Hz so 4 secs window is 40 datapoints
config.Speed.velocityCutoff = 0.2;                      % velocity cutoff to select only the walking part of the reconstructed velocity
config.Speed.timeOffsetForDetectedTemporalWindow = 0.5; % time in seconds that will push earlier/ the detected rising edge

config.ModelName = "LIModel";
config.NumParams = 3;

resultfolder = savefolder+config.ModelName+"/";
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end
config.ResultFolder = resultfolder;

config.TrialFilter = 0; %merge all conditions
config.UseGlobalSearch = true;

%% calculating tracking path and transoform data
disp('%%%%%%%%%%%%%%% Calculating Tracking Path - YOUNG ... %%%%%%%%%%%%%%%');
% threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
config.Speed.tresholdForBadParticipantL1Recontruction = 1.55;
YoungControls   = CalculateTrackingPath(YoungControls, config);
%transform data
YoungControls = TransformPaths(YoungControls);

%% 
YoungResults = PerformGroupFit(YoungControls, config);



