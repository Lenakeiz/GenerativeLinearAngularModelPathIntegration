%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default'); %for code reproducibility

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
savefolder = pwd + "/Output/";

%% setting the configuration
config.Speed.alpha                                  = 0.9;                 %Paramanter for running speed calculation
config.Speed.timeOffsetAfterFlagReach               = 1.5;                 %Time to track after flag reached in seconds 
config.Speed.smoothWindow                           = 10;                  % tracking rate should be 10Hz so 4 secs window is 40 datapoints
config.Speed.velocityCutoff                         = 0.2;                 % velocity cutoff to select only the walking part of the reconstructed velocity
config.Speed.timeOffsetForDetectedTemporalWindow    = 0.4;                 % time in seconds that will push earlier/ the detected rising edge
config.UseGlobalSearch                              = true;
config.TrackedInboundAngularDeltaT                  = 1;
config.includeStand                                 = false;
config.useOoBTrial                                  = false;
config.useweber                                     = false; %only true when use weber law in simple generative models
%% Model fitting

config.ModelName        = "ConstSpeedModelDistAngleRGb";
config.ParamName        = ["beta", "sigma", "g3", "b", "nu"];
config.NumParams        = length(config.ParamName);

%% Model fitting for HealthyOld
config.Speed.tresholdForBadParticipantL1Recontruction = 2.0; 
HealthyControls = TransformPaths(HealthyControls);
HealthyControls   = CalculateTrackingPath(HealthyControls, config);
ManuallyScoringHealthyOld;

%% Model fitting for MCIPos
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCIPos = TransformPaths(MCIPos);%transform data
MCIPos   = CalculateTrackingPath(MCIPos, config);
ManuallyScoringMCIPos;
[AllMCIPosParams, AllMCIPosX, AllMCIPosDX, AllMCIPosTheta, AllMCIPosIC] = getResultsAllConditions(MCIPos, config);

%% Model fitting for MCINeg
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCINeg = TransformPaths(MCINeg);%transform data
MCINeg   = CalculateTrackingPath(MCINeg, config);
ManuallyScoringMCINeg;
%[AllMCINegParams,  ~, ~, ~, ~] = getResultsAllConditions(MCINeg, config);

%% Model fitting for MCIUnk
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCIUnk = TransformPaths(Unknown);%transform data
MCIUnk   = CalculateTrackingPath(MCIUnk, config);
ManuallyScoringMCIUnk;
%[AllMCINegParams,  ~, ~, ~, ~] = getResultsAllConditions(MCINeg, config);
