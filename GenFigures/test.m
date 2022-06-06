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
config.TrackedInboundAngularDeltaT                  = 5;
config.UseGlobalSearch = true;
resultfolder = savefolder+"PaperFigs/FigVisualization";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% Model fitting for YoungControl data
%% calculating tracking path and transoform data
config.Speed.tresholdForBadParticipantL1Recontruction = 1.55;   % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
YoungControls   = TransformPaths(YoungControls); %transform data
YoungControls   = CalculateTrackingPath(YoungControls, config);
%VisualizeRealtimeTrackingData(YoungControls,1, 1, 0.0001, 'cutconethree',true);

%% Model fitting for HealthyOld data
config.Speed.tresholdForBadParticipantL1Recontruction = 2.0; 
HealthyControls = TransformPaths(HealthyControls);%transform data
HealthyControls = CalculateTrackingPath(HealthyControls, config);

%VisualizeRealtimeTrackingData(HealthyControls,1, 1, 0.0001, 'cutconethree',true);
%% Model fitting for MCIPos
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCIPos          = TransformPaths(MCIPos);%transform data
MCIPos          = CalculateTrackingPath(MCIPos, config);

%% Model fitting for MCINeg
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCINeg          = TransformPaths(MCINeg);%transform data
MCINeg          = CalculateTrackingPath(MCINeg, config);
%VisualizeRealtimeTrackingData(MCINeg,1, 1, 0.0001, 'cutconethree',true);

%% Model fitting for MCIUnk
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
Unknown          = TransformPaths(Unknown);%transform data
Unknown         = CalculateTrackingPath(Unknown, config);

%%
VisualizeRealtimeTrackingData(MCIUnk,1, 1, 0.0001, 'cutconethree',true);
%% Setting colors for using in plots
ColorPattern; 