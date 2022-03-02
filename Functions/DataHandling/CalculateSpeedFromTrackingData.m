%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default');

disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
    
%%

YoungControls   = CalculateTrackingPath(YoungControls);
HealthyControls = CalculateTrackingPath(HealthyControls);
Unknown         = CalculateTrackingPath(Unknown);
MCIPos          = CalculateTrackingPath(MCIPos);
MCINeg          = CalculateTrackingPath(MCINeg);