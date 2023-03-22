%% Vector Addiction Modelling (VAM) Model
% Andrea Castegnaro, UCL, 2022 uceeaca@ucl.ac.uk
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