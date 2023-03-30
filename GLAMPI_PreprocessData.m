%% Generative Linear-Angular Model of Path Integration (GLAMPI) Preprocessing stage
% Andrea Castegnaro, UCL, 2022 uceeaca@ucl.ac.uk
% This script takes the loaded data and perform the preprocessing stage as
% described in the online methods for all groups of participants.

%% Preprocessing Young
disp("%%%%%%%%%%%%%%% Preprocessing young controls %%%%%%%%%%%%%%%");
% In early version of the task the algorithm for creating the triangle was
% not refined to maximize the walking length, resulting in some triangles
% that were much smaller than the others. We want to flag these trials and
% exclude participants that used an early version of the task.
config.Speed.tresholdForBadParticipantL1Recontruction   = 1.55;    
YoungControls = TransformPaths(YoungControls); %transform data
YoungControls   = CalculateTrackingPath(YoungControls, config);
ManuallyScoringYoung;
%% Preprocessing Elderly
disp("%%%%%%%%%%%%%%% Preprocessing healthy controls %%%%%%%%%%%%%%%");
% In early version of the task the algorithm for creating the triangle was
% not refined to maximize the walking length, resulting in some triangles
% that were much smaller than the others. We want to flag these trials and
% exclude participants that used an early version of the task.
config.Speed.tresholdForBadParticipantL1Recontruction   = 2.0;    % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
HealthyControls = TransformPaths(HealthyControls);%transform data
HealthyControls   = CalculateTrackingPath(HealthyControls, config);
ManuallyScoringHealthyOld;
%% Preprocessing MCI positive
disp("%%%%%%%%%%%%%%% Starting fit for mci positive %%%%%%%%%%%%%%%");
% All of the mci have been tested with the final version of the algorithm
% creation.
config.Speed.tresholdForBadParticipantL1Recontruction   = 0.0;
MCIPos = TransformPaths(MCIPos);%transform data
MCIPos   = CalculateTrackingPath(MCIPos, config);
ManuallyScoringMCIPos;
%% Preprocessing MCI negative
disp("%%%%%%%%%%%%%%% Preprocessing mci negative %%%%%%%%%%%%%%%");
config.Speed.tresholdForBadParticipantL1Recontruction   = 0.0; 
MCINeg = TransformPaths(MCINeg);%transform data
MCINeg   = CalculateTrackingPath(MCINeg, config);
ManuallyScoringMCINeg;
%% Preprocessing MCI unknown
disp("%%%%%%%%%%%%%%% Preprocessing mci unknown %%%%%%%%%%%%%%%");
config.Speed.tresholdForBadParticipantL1Recontruction   = 0.0; 
MCIUnk = TransformPaths(Unknown);%transform data
MCIUnk   = CalculateTrackingPath(MCIUnk, config);
ManuallyScoringMCIUnk;
clear Unknown % removing duplicate

disp("%%%%%%%%%%%%%%% Preprocessing complete %%%%%%%%%%%%%%%");
