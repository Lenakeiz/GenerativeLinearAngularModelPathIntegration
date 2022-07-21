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

if (~isfield(config,'ModelName') | ~isfield(config,'ParamName') | ~isfield(config,'NumParams'))
    ME = MException("MyComponent:noSuchVariable","Model variables must be specified before running the model: set ModelName, ParamName and NumParams in config struct");
    throw(ME);
end

%% Model fitting for Young Controls
disp("%%%%%%%%%%%%%%% Starting fit for young controls %%%%%%%%%%%%%%%");
YoungControls.Results = getResultsAllConditions(YoungControls, config);

%% Model fitting for Healthy Controls
disp("%%%%%%%%%%%%%%% Starting fit for healthy controls %%%%%%%%%%%%%%%");
HealthyControls.Results = getResultsAllConditions(HealthyControls, config);

%% Model fitting for MCI Pos
disp("%%%%%%%%%%%%%%% Starting fit for mci pos %%%%%%%%%%%%%%%");
MCIPos.Results = getResultsAllConditions(MCIPos, config);

%% Model fitting for MCI Neg
disp("%%%%%%%%%%%%%%% Starting fit for mci neg %%%%%%%%%%%%%%%");
MCINeg.Results = getResultsAllConditions(MCINeg, config);

%% Model fitting for MCI Unk
disp("%%%%%%%%%%%%%%% Starting fit for mci unknown %%%%%%%%%%%%%%%");
MCIUnk.Results = getResultsAllConditions(MCIUnk, config);

disp("%%%%%%%%%%%%%%% Model fitting complete %%%%%%%%%%%%%%%");