%% Prepare a simple config file for the model to run
%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default'); %for code reproducibility

%% Loading data
disp('%%%%%%%%%%%%%%% Data loading ... %%%%%%%%%%%%%%%');
load('Data/HowettBrain2019_Dataset.mat');

%% Loading MRI data
disp('%%%%%%%%%%%%%%% Adding MRI data ... %%%%%%%%%%%%%%%');
LoadMRIData;

%% Loading NeuroPsychological data
disp('%%%%%%%%%%%%%%% Adding NeuroPsychological data ... %%%%%%%%%%%%%%%');
LoadNeuroPsychData;

%% setting the configuration
config.Speed.alpha                                      = 0.9;    % Paramanter for running speed calculation
config.Speed.timeOffsetAfterFlagReach                   = 1.5;    % Time to track after flag reached in seconds 
config.Speed.smoothWindow                               = 10;     % tracking rate should be 10Hz so 4 secs window is 40 datapoints
config.Speed.velocityCutoff                             = 0.2;    % velocity cutoff to select only the walking part of the reconstructed velocity
config.Speed.timeOffsetForDetectedTemporalWindow        = 0.4;    % time in seconds that will push earlier/ the detected rising edge
config.UseGlobalSearch                                  = true;
config.TrackedInboundAngularDeltaT                      = 1;
config.includeStand                                     = false;
config.useweber                                         = false;  % only true when use weber law in simple generative models
config.useOoBtrials                                     = true;
config.Speed.tresholdForBadParticipantL1Recontruction   = 0.0;
config.useTrialFilter                                   = true;   % when true the model will be fitted for each of the task conditions separately. If false it will discard the

disp('%%%%%%%%%%%%%%% Data loading complete... %%%%%%%%%%%%%%%');