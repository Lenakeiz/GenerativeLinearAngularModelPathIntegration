%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default');

disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');

%%

config.Speed.alpha = 0.9; %Paramanter for running speed calculation
config.Speed.TOffsetAfterFlagReach = 2; %Time to track after flag reached in seconds 
config.Speed.smoothWindow = 40; % tracking rate should be 10Hz so 4 secs window

YoungControls   = CalculateTrackingPath(YoungControls, config);
%% 
HealthyControls = CalculateTrackingPath(HealthyControls, config);
Unknown         = CalculateTrackingPath(Unknown, config);
MCIPos          = CalculateTrackingPath(MCIPos, config);
MCINeg          = CalculateTrackingPath(MCINeg, config);

%% Plotting cumulative plots

pId = 5;

nTrials = size(YoungControls.TrackedL1{1,pId},1);

figure;

subplot("Position", [0.05 0.05 0.4 0.8]);
hold on;

for tr = 1:nTrials
    
    maxL1 = YoungControls.TrackedL1{1,pId}{tr,1}.Smoothed_Vel_proj;
    idxMaxL1 = find(YoungControls.TrackedL1{1,pId}{tr,1}.Smoothed_Vel_proj == maxL1);
    
    tempTime = YoungControls.TrackedL1{1,pId}{tr,1}.Time - YoungControls.TrackedL1{1,pId}{tr,1}.Time(idxMaxL1);

    plot(tempTime, YoungControls.TrackedL1{1,pId}{tr,1}.Smoothed_Vel_proj,...Color=[config.color_scheme_npg(9,:) 1],...
        LineStyle="-",...
        LineWidth=1.3);

end

hold off;

%subplot("Position", [0.55 0.05 0.4 0.8]);