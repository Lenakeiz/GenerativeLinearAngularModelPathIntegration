%% Main entry point for the Vector Addiction Modelling (VAM) project
% This script loads the data the from the Howett 2019., Brain paper. It
% prepare the variables and runs the modelling described in the current
% paper. Each figure will have a separate script that can be run after this
% script has been correctly executed. For additional dependecies please
% refer to the README file.

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/HowettBrain2019_Dataset.mat');

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
config.useOoBtrials = true;

%% Model fitting for Young Controls
disp("%%%%%%%%%%%%%%% Starting fit for young controls %%%%%%%%%%%%%%%");
config.Speed.tresholdForBadParticipantL1Recontruction = 1.55;  
YoungControls = TransformPaths(YoungControls);%transform data
YoungControls   = CalculateTrackingPath(YoungControls, config);
ManuallyScoringYoung;
YoungControls.Results = getResultsAllConditions(YoungControls, config);

%% Model fitting for Healthy Controls
disp("%%%%%%%%%%%%%%% Starting fit for healthy controls %%%%%%%%%%%%%%%");
config.Speed.tresholdForBadParticipantL1Recontruction = 2.0;
HealthyControls = TransformPaths(HealthyControls);%transform data
HealthyControls   = CalculateTrackingPath(HealthyControls, config);
ManuallyScoringHealthyOld;
HealthyControls.Results = getResultsAllConditions(HealthyControls, config);

%% Model fitting for MCI Pos
disp("%%%%%%%%%%%%%%% Starting fit for mci positive %%%%%%%%%%%%%%%");
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCIPos = TransformPaths(MCIPos);%transform data
MCIPos   = CalculateTrackingPath(MCIPos, config);
ManuallyScoringMCIPos;
MCIPos.Results = getResultsAllConditions(MCIPos, config);

%% Model fitting for MCI Neg
disp("%%%%%%%%%%%%%%% Starting fit for mci negative %%%%%%%%%%%%%%%");
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCINeg = TransformPaths(MCINeg);%transform data
MCINeg   = CalculateTrackingPath(MCINeg, config);
ManuallyScoringMCINeg;
MCINeg.Results = getResultsAllConditions(MCINeg, config);

%% Model fitting for MCI Unk
disp("%%%%%%%%%%%%%%% Starting fit for mci unknown %%%%%%%%%%%%%%%");
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCIUnk = TransformPaths(Unknown);%transform data
MCIUnk   = CalculateTrackingPath(MCIUnk, config);
ManuallyScoringMCIUnk;
MCIUnk.Results = getResultsAllConditions(MCIUnk, config);