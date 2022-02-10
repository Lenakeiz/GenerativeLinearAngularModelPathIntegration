%% FittingData 
% This script will load the data, rotate the paths and call the solver to find the parameters 
% that minimize the mean distance to generated x3' which matchs the physical leg 3

%% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
%load('Data/AllDataErrors2018_V2.mat');
%V3 including the missing participants, which should always being loaded now
load('Data/AllDataErrors2018_V3.mat');

savefolder = "C:/Users/Zilong/Desktop/path integration model/Andrea's matlab code/PathIntegrationDataEstimationModelling/Output/";


%% choose the model that we want to optimize

%choose the model that we want to optimize
%the base model with G1, G2, G3, g2, g3 k3, sigma, nu. #params=8
%Model_Name = "base_add_G3"; numParams=8;  

%the base model with G1, G2, G3=1, g2, g3, k3, sigma, nu. #params=7
Model_Name = "base"; numParams=7;   

%consider k3=0 in the excution turning error of leg 3, i.e., G1, G2, G3=1, g2, g3, k3=0, sigma, nu. #params=6
%Model_Name = "set_k3_0"; numParams=6;   

%no consider of execution error i.e., G1, G2, G3=1, g2, g3=1, k3=0, sigma, nu. #params=5
%Model_Name = "set_g3_1_k3_0"; numParams=5; 

%no consider of the rotation gain factor in leg2, i.e., G1, G2, G3=1, g2=1, g3, k3, sigma, nu. #params=6
%Model_Name = "set_g2_1"; numParams=6;  

%no consider of the length gain factor in leg2, i.e., G1, G2=1, G3=1, g2, g3, k3, sigma, nu. #params=6
%Model_Name = "set_bG2_1"; numParams=6;   

%no consider of the length gain factor in leg1, i.e., G1=1, G2, G3=1, g2, g3, k3, sigma, nu. #params=6
%Model_Name = "set_bG1_1"; numParams=6;   

%same length gain factor in leg1 & 2, i.e., G1=G2, g2, g3, k3, sigma, nu. #params=6
%Model_Name = "set_same_bG"; numParams=6; 

%same rotation gain factor in leg2 & 3, i.e., G1, G2, g2=g3, k3=0, sigma, nu.#params=5
%Model_Name = "set_same_sg"; numParams=5; 

%same length gain, same rotation gain, i.e., G1=G2, g2=g3, k3=0, sigma, nu.#params=4
%Model_Name = "set_same_Gg"; numParams=4; 


%% initial a empty cell for storing all configuration (for passing varibles to functions easliy) 
config.ModelName = Model_Name;
config.NumParams = numParams;

config.USEOoBTrials = true;
if config.USEOoBTrials==true
    resultfolder = savefolder+Model_Name+"/OoB/";
else
    resultfolder = savefolder+Model_Name+"/NoOoB/";
end
config.ResultFolder = resultfolder;

%create storing folder for store data if not exist
config.TrialFilter = 1;
config.TrialFolder = resultfolder+"TrailType"+num2str(config.TrialFilter)+"/";
if ~exist(config.TrialFolder, 'dir')
   mkdir(config.TrialFolder);
end

%plot
config.ifReturnAngleDistribution = true;
config.ifPhysicalTrajectory = true;
config.ifMentalTrajectory = true;
config.ifEndpointDistribution = true;

%% storing the bayesian inference criterion in a cell
IC_ALL = cell(1,5);

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Performing data fitting on Young participants
%load the data
if config.USEOoBTrials==true
    YoungControls = TransformPathsOoB(YoungControls);
else
    YoungControls = TransformPaths(YoungControls);
end
%do fitting
disp('%%%%%%%%%%%%%%% PERFORMING FIT YOUNG CONTROLS %%%%%%%%%%%%%%%');
YoungResults = PerformGroupFit(YoungControls, config);
IC_ALL{1} = YoungResults.IC;
if config.ifReturnAngleDistribution==true
    % plot return angle turn distribution including alpha, theta3prime etc
    PlotReturnAngleDistribution(YoungResults, 'Young', config);
end
if config.ifPhysicalTrajectory==true
    % plot physical trajectories
    PlotPhysicalTrajectory(YoungResults, 'Young', config);
end
if config.ifMentalTrajectory==true
    % plot mental trajectories
    PlotMentalTrajectory(YoungResults, 'Young', config);
end
if config.ifEndpointDistribution==true
    % plot end points distribution
    PlotEndpointDistribution(YoungResults, 'Young', config);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Performing data fitting on HealthyOld participants
%load the data
if config.USEOoBTrials==true
    HealthyOld = TransformPathsOoB(HealthyControls);
else
    HealthyOld = TransformPaths(HealthyControls);
end
%do fitting
disp('%%%%%%%%%%%%%%% PERFORMING FIT OLDER CONTROLS %%%%%%%%%%%%%%%');
HealthyOldResults = PerformGroupFit(HealthyOld, config);
IC_ALL{2} = HealthyOldResults.IC; 
if config.ifReturnAngleDistribution==true
    % plot return angle turn distribution including alpha, theta3prime etc
    PlotReturnAngleDistribution(HealthyOldResults, 'HealthyOld', config);
end
if config.ifPhysicalTrajectory==true
    % plot physical trajectories
    PlotPhysicalTrajectory(HealthyOldResults, 'HealthyOld', config);
end
if config.ifMentalTrajectory==true
    % plot mental trajectories
    PlotMentalTrajectory(HealthyOldResults, 'HealthyOld', config);
end
if config.ifEndpointDistribution==true
    % plot end points distribution
    PlotEndpointDistribution(HealthyOldResults, 'HealthyOld', config);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Performing data fitting on MCIPos participants
%load the data
if config.USEOoBTrials==true
    MCIPos = TransformPathsOoB(MCIPos);
else
    MCIPos = TransformPaths(MCIPos);
end
%do fitting
disp('%%%%%%%%%%%%%%% PERFORMING FIT MCIPos %%%%%%%%%%%%%%%');
MCIPosResults = PerformGroupFit(MCIPos, config);
IC_ALL{3} = MCIPosResults.IC; 

if config.ifReturnAngleDistribution==true
    % plot return angle turn distribution including alpha, theta3prime etc
    PlotReturnAngleDistribution(MCIPosResults, 'MCIPos', config);
end
if config.ifPhysicalTrajectory==true
    % plot physical trajectories
    PlotPhysicalTrajectory(MCIPosResults, 'MCIPos', config);
end
if config.ifMentalTrajectory==true
    % plot mental trajectories
    PlotMentalTrajectory(MCIPosResults, 'MCIPos', config);
end
if config.ifEndpointDistribution==true
    % plot end points distribution
    PlotEndpointDistribution(MCIPosResults, 'MCIPos', config);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Performing data fitting on MCINeg participants
%load the data
if config.USEOoBTrials==true
    MCINeg = TransformPathsOoB(MCINeg);
else
    MCINeg = TransformPaths(MCINeg);
end
%do fitting
disp('%%%%%%%%%%%%%%% PERFORMING FIT MCINeg %%%%%%%%%%%%%%%');
MCINegResults = PerformGroupFit(MCINeg, config);
IC_ALL{4} = MCINegResults.IC; 

if config.ifReturnAngleDistribution==true
    % plot return angle turn distribution including alpha, theta3prime etc
    PlotReturnAngleDistribution(MCINegResults, 'MCINeg', config);
end
if config.ifPhysicalTrajectory==true
    % plot physical trajectories
    PlotPhysicalTrajectory(MCINegResults, 'MCINeg', config);
end
if config.ifMentalTrajectory==true
    % plot mental trajectories
    PlotMentalTrajectory(MCINegResults, 'MCINeg', config);
end
if config.ifEndpointDistribution==true
    % plot end points distribution
    PlotEndpointDistribution(MCINegResults, 'MCINeg', config);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Performing data fitting on MCIUnk participants
%load the data
if config.USEOoBTrials==true
    MCIUnk = TransformPathsOoB(Unknown);
else
    MCIUnk = TransformPaths(Unknown);
end
%do fitting
disp('%%%%%%%%%%%%%%% PERFORMING FIT MCIUnk %%%%%%%%%%%%%%%');
MCIUnkResults = PerformGroupFit(MCIUnk, config);
IC_ALL{5} = MCIUnkResults.IC;
if config.ifReturnAngleDistribution==true
    % plot return angle turn distribution including alpha, theta3prime etc
    PlotReturnAngleDistribution(MCIUnkResults, 'MCIUnk', config);
end
if config.ifPhysicalTrajectory==true
    % plot physical trajectories
    PlotPhysicalTrajectory(MCIUnkResults, 'MCIUnk', config);
end
if config.ifMentalTrajectory==true
    % plot mental trajectories
    PlotMentalTrajectory(MCIUnkResults, 'MCIUnk', config);
end
if config.ifEndpointDistribution==true
    % plot end points distribution
    PlotEndpointDistribution(MCIUnkResults, 'MCIUnk', config);
end
save(resultfolder+"IC_"+Model_Name+".mat", 'IC_ALL');

%% performing one-way anova
OnewayAnovaOnParams(YoungResults.estimatedParams, ...
                    HealthyOldResults.estimatedParams, ...
                    MCIPosResults.estimatedParams, ...
                    MCINegResults.estimatedParams, ...
                    MCIUnkResults.estimatedParams, config);

%% error bar plot of estimated parameters across participants
PlotErrorBar(YoungResults.estimatedParams, ...
            HealthyOldResults.estimatedParams, ...
            MCIPosResults.estimatedParams, ...
            MCINegResults.estimatedParams, ...
            MCIUnkResults.estimatedParams, config);