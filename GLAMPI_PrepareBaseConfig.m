%% Generative Linear-Angular Model of Path Integration (GLAMPI) Base Configurator
% Andrea Castegnaro, UCL, 2022 uceeaca@ucl.ac.uk
% This script loads the Howett et al., 2019 behavioural and MRI data. The
% dataset contains also the healthy young participants collected elsewhere.
% It also set global properties that are used for the preprocecssing
% stage and for fitting the model. Information for the config parameters
% can be found within the script.

%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default');

%% Loading data
disp('%%%%%%%%%%%%%%% Data loading ... %%%%%%%%%%%%%%%');
load('Data/HowettBrain2019_Dataset.mat');

%% setting the configuration
config.Speed.alpha                                      = 0.9;    % Paramanter for running speed calculation
config.Speed.timeOffsetAfterFlagReach                   = 0.5;    % Time to track after cone reached in seconds 
config.Speed.smoothWindow                               = 20;     % tracking rate should be 10Hz so 4 secs window is 40 datapoints
config.Speed.velocityCutoff                             = 0.2;    % velocity cutoff to select only the walking part of the reconstructed velocity
config.Speed.timeOffsetForDetectedTemporalWindow        = 0.2;    % time in seconds that will push earlier/ the detected rising edge
config.TrackedInboundAngularDeltaT                      = 1;      % delta time step to integrate the angular information from the tracking data
config.includeStand                                     = false;  % set to true if the model should also include for duration of time spent on each segment the time calculate from tracking data while standing at each cone
config.useOoBtrials                                     = true;   % Are we using out of bound trials in the model?
config.useTrialFilter                                   = true;   % when true the model will be fitted for each of the task conditions separately. If false it will discard the

disp('%%%%%%%%%%%%%%% Data loading complete... %%%%%%%%%%%%%%%');