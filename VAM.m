%% Main entry point for the Vector Addiction Modelling (VAM) project
% This script loads the data the from the Howett 2019., Brain paper. It
% prepare the variables and runs the modelling described in the current
% paper. Each figure will have a separate script that can be run after this
% script has been correctly executed. For additional dependecies please
% refer to the README file. 
% This script requires a config variable set before running. An example of
% the config variable can be created with VAM_PrepareBaseConfig

if (exist('config','var') == 0)
    ME = MException("MyComponent:noSuchVariable","Config variable not found, please run VAM_PrepareBaseConfig to generate a base one");
    throw(ME);
end

%% Model fitting for Young Controls
config.Speed.tresholdForBadParticipantL1Recontruction   = 1.55;    % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
disp("%%%%%%%%%%%%%%% Starting fit for young controls %%%%%%%%%%%%%%%");
YoungControls = TransformPaths(YoungControls);%transform data
YoungControls   = CalculateTrackingPath(YoungControls, config);
ManuallyScoringYoung;
YoungControls.Results = getResultsAllConditions(YoungControls, config);

%% Model fitting for Healthy Controls
config.Speed.tresholdForBadParticipantL1Recontruction   = 2.0;    % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
disp("%%%%%%%%%%%%%%% Starting fit for healthy controls %%%%%%%%%%%%%%%");
HealthyControls = TransformPaths(HealthyControls);%transform data
HealthyControls   = CalculateTrackingPath(HealthyControls, config);
ManuallyScoringHealthyOld;
HealthyControls.Results = getResultsAllConditions(HealthyControls, config);

%% Model fitting for MCI Pos
config.Speed.tresholdForBadParticipantL1Recontruction   = 0.0;    % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
disp("%%%%%%%%%%%%%%% Starting fit for mci positive %%%%%%%%%%%%%%%");
MCIPos = TransformPaths(MCIPos);%transform data
MCIPos   = CalculateTrackingPath(MCIPos, config);
ManuallyScoringMCIPos;
MCIPos.Results = getResultsAllConditions(MCIPos, config);

%% Model fitting for MCI Neg
disp("%%%%%%%%%%%%%%% Starting fit for mci negative %%%%%%%%%%%%%%%");
MCINeg = TransformPaths(MCINeg);%transform data
MCINeg   = CalculateTrackingPath(MCINeg, config);
ManuallyScoringMCINeg;
MCINeg.Results = getResultsAllConditions(MCINeg, config);

%% Model fitting for MCI Unk
disp("%%%%%%%%%%%%%%% Starting fit for mci unknown %%%%%%%%%%%%%%%");
MCIUnk = TransformPaths(Unknown);%transform data
MCIUnk   = CalculateTrackingPath(MCIUnk, config);
ManuallyScoringMCIUnk;
MCIUnk.Results = getResultsAllConditions(MCIUnk, config);