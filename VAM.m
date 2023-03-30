%% Generative Linear-Angular Model of Path Integration (GLAMPI)
% Andrea Castegnaro, UCL, 2022 uceeaca@ucl.ac.uk
% Entry point for data analysis and modelling. This script should be
% executed after the preprocessing step and it is a wrapper script. For
% details on how to use it please refer to any of the script inside the
% Main folder

if (exist('config','var') == 0)
    ME = MException("MyComponent:noSuchVariable","Config variable not found, please run VAM_PrepareBaseConfig to generate a base one");
    throw(ME);
end

if (~isfield(config,'ModelName') | ~isfield(config,'ParamName') | ~isfield(config,'NumParams'))
    ME = MException("MyComponent:noSuchVariable","Model variables must be specified before running the model: set ModelName, ParamName and NumParams in config struct");
    throw(ME);
end

%% Using a dummy variable "Model Selection" to check if we need to run the model only the Healthy elderly
% Model fitting for Young Controls
if (isfield(config, 'ModelSelection') == 0)
    disp("%%%%%%%%%%%%%%% Starting fit for young controls %%%%%%%%%%%%%%%");
    YoungControls.Results = getResultsAllConditions(YoungControls, config);
end
% Model fitting for Healthy Controls
disp("%%%%%%%%%%%%%%% Starting fit for healthy controls %%%%%%%%%%%%%%%");
HealthyControls.Results = getResultsAllConditions(HealthyControls, config);

% Model fitting for MCI Pos
if (isfield(config, 'ModelSelection') == 0)
    disp("%%%%%%%%%%%%%%% Starting fit for mci pos %%%%%%%%%%%%%%%");
    MCIPos.Results = getResultsAllConditions(MCIPos, config);
end
% Model fitting for MCI Neg
if (isfield(config, 'ModelSelection') == 0)
    disp("%%%%%%%%%%%%%%% Starting fit for mci neg %%%%%%%%%%%%%%%");
    MCINeg.Results = getResultsAllConditions(MCINeg, config);
end
% Model fitting for MCI Unk
if (isfield(config, 'ModelSelection') == 0)
    disp("%%%%%%%%%%%%%%% Starting fit for mci unknown %%%%%%%%%%%%%%%");
    MCIUnk.Results = getResultsAllConditions(MCIUnk, config);
end

disp("%%%%%%%%%%%%%%% Model fitting complete %%%%%%%%%%%%%%%");