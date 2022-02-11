%% FittingData 
% This script will lo%% FittingData 
% This script will load the data, rotate the paths and call the solver to find the parameters 
% that minimize the mean distance to generated x3' which matchs the physical leg 3

%% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
%load('Data/AllDataErrors2018_V2.mat');
%V3 including the missing participants, which should always being loaded now
load('Data/AllDataErrors2018_V3.mat');
savefolder = "C:/Users/Zilong/Desktop/path integration model/VectorAdditionModel/Output/";

%% choose the model 

%the full model with G1, G2, G3, g2, g3 k3, sigma, nu. #params=8
%Model_Name = "base_add_G3"; numParams=8;  

%the base model with G1, G2, G3=1, g2, g3, k3, sigma, nu. #params=7
%Model_Name = "base"; numParams=7;   

%consider k3=0 in the excution turning error of leg 3, i.e., G1, G2, G3=1, g2, g3, k3=0, sigma, nu. #params=6
%Model_Name = "set_k3_0"; numParams=6;   

%no consider of execution error i.e., G1, G2, G3=1, g2, g3=1, k3=0, sigma, nu. #params=5
%Model_Name = "set_g3_1_k3_0"; numParams=5; 

%no consider of the rotation gain factor in leg2, i.e., G1, G2, G3=1, g2=1, g3, k3, sigma, nu. #params=6
Model_Name = "BaseModel"; numParams=5;  

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
    resultfolder = savefolder+"AllGroups/"+Model_Name+"/OoB/";
else
    resultfolder = savefolder+"AllGroups/"+Model_Name+"/NoOoB/";
end
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end
config.ResultFolder = resultfolder;
%
config.flagPlotCorrelation = false;
config.flagPlotReturnAngleDistribution = false;

%for storing all IC information
%create storing folder for trajectory if not exist
icfolder = resultfolder+"IC/";
if ~exist(icfolder, 'dir')
   mkdir(icfolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over three different conditions for Young participants
%load the data
if config.USEOoBTrials==true
    YoungControls = TransformPathsOoB(YoungControls);
else
    YoungControls = TransformPaths(YoungControls);
end
%get fitting results for all conditions
[AllYoungParams, AllYoungX, AllYoungDX, AllYoungTheta, AllYoungflagOoB, AllYoungIC] = getResultsAllConditions(YoungControls, config);
%storing IC
save(icfolder+"ICYoung.mat", 'AllYoungIC');
%%
if config.flagPlotCorrelation == true
    %plot correlation between turn angle and k3
    AllYoungTurnAngleDiff = PlotCorrelationBetweenReturnAngleAndK3(AllYoungParams, AllYoungDX, AllYoungTheta, 'Young', config);
end
%%
if config.flagPlotReturnAngleDistribution
    %plot last turn angle distribution
    PlotReturnAngleDistribution(AllYoungParams, AllYoungDX, AllYoungTheta, AllYoungflagOoB, 'Young', config);
    PlotReturnAngleDistributionInVerticalScatter(AllYoungParams, AllYoungDX, AllYoungTheta, 'Young', config);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over three different conditions for HeathyOld participants
%load the data
if config.USEOoBTrials==true
    HealthyOld = TransformPathsOoB(HealthyControls);
else
    HealthyOld = TransformPaths(HealthyControls);
end
%get fitting results for all conditions
[AllHealthyOldParams, AllHealthyOldX, AllHealthyOldDX, AllHealthyOldTheta, AllHealthyOldflagOoB, AllHealthyOldIC] = getResultsAllConditions(HealthyOld, config);
%storing IC
save(icfolder+"ICHealthyOld.mat", 'AllHealthyOldIC');
if config.flagPlotCorrelation == true
    %plot correlation between turn angle and k3
    AllHealthyOldTurnAngleDiff = PlotCorrelationBetweenReturnAngleAndK3(AllHealthyOldParams, AllHealthyOldDX, AllHealthyOldTheta, 'HealthyOld', config);
end
%%
if config.flagPlotReturnAngleDistribution
    %plot last turn angle distribution
    PlotReturnAngleDistribution(AllHealthyOldParams, AllHealthyOldDX, AllHealthyOldTheta, AllHealthyOldflagOoB, 'HealthyOld', config);
    PlotReturnAngleDistributionInVerticalScatter(AllHealthyOldParams, AllHealthyOldDX, AllHealthyOldTheta, 'HealthyOld', config);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over three different conditions for MCIPos participants
%load the data
if config.USEOoBTrials==true
    MCIPos = TransformPathsOoB(MCIPos);
else
    MCIPos = TransformPaths(MCIPos);
end
%get fitting results for all conditions
[AllMCIPosParams, AllMCIPosX, AllMCIPosDX, AllMCIPosTheta, AllMCIPosflagOoB, AllMCIPosIC] = getResultsAllConditions(MCIPos, config);
%storing IC
save(icfolder+"ICMCIPos.mat", 'AllMCIPosIC');
if config.flagPlotCorrelation == true
    %plot correlation between turn angle and k3
    AllMCIPosTurnAngleDiff = PlotCorrelationBetweenReturnAngleAndK3(AllMCIPosParams, AllMCIPosDX, AllMCIPosTheta, 'MCIPos', config);
end
%%
if config.flagPlotReturnAngleDistribution
    %plot last turn angle distribution
    PlotReturnAngleDistribution(AllMCIPosParams, AllMCIPosDX, AllMCIPosTheta, AllMCIPosflagOoB, 'MCIPos', config);
    PlotReturnAngleDistributionInVerticalScatter(AllMCIPosParams, AllMCIPosDX, AllMCIPosTheta, 'MCIPos', config);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over three different conditions for MCINeg participants
%load the data
if config.USEOoBTrials==true
    MCINeg = TransformPathsOoB(MCINeg);
else
    MCINeg = TransformPaths(MCINeg);
end
%get fitting results for all conditions
[AllMCINegParams, AllMCINegX, AllMCINegDX, AllMCINegTheta, AllMCINegflagOoB, AllMCINegIC] = getResultsAllConditions(MCINeg, config);
%storing IC
save(icfolder+"ICMCINeg.mat", 'AllMCINegIC');
if config.flagPlotCorrelation == true
    %plot correlation between turn angle and k3
    AllMCINegTurnAngleDiff = PlotCorrelationBetweenReturnAngleAndK3(AllMCINegParams, AllMCINegDX, AllMCINegTheta, 'MCINeg', config);
end
%%
if config.flagPlotReturnAngleDistribution
    %plot last turn angle distribution
    PlotReturnAngleDistribution(AllMCINegParams, AllMCINegDX, AllMCINegTheta, AllMCINegflagOoB, 'MCINeg', config);
    PlotReturnAngleDistributionInVerticalScatter(AllMCINegParams, AllMCINegDX, AllMCINegTheta, 'MCINeg', config);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over three different conditions for MCIUnk participants
%load the data
if config.USEOoBTrials==true
    MCIUnk= TransformPathsOoB(Unknown);
else
    MCIUnk = TransformPaths(Unknown);
end
%get fitting results for all conditions
[AllMCIUnkParams, AllMCIUnkX, AllMCIUnkDX, AllMCIUnkTheta, AllMCIUnkflagOoB, AllMCIUnkIC] = getResultsAllConditions(MCIUnk, config);
%storing IC
save(icfolder+"ICMCIUnk.mat", 'AllMCIUnkIC');
if config.flagPlotCorrelation == true
    %plot correlation between turn angle and k3
    AllMCIUnkTurnAngleDiff = PlotCorrelationBetweenReturnAngleAndK3(AllMCIUnkParams, AllMCIUnkDX, AllMCIUnkTheta, 'MCIUnk', config);
end
%%
if config.flagPlotReturnAngleDistribution
    %plot last turn angle distribution
    PlotReturnAngleDistribution(AllMCIUnkParams, AllMCIUnkDX, AllMCIUnkTheta, AllMCIUnkflagOoB, 'MCIUnk', config);
    PlotReturnAngleDistributionInVerticalScatter(AllMCIUnkParams, AllMCIUnkDX, AllMCIUnkTheta, 'MCIUnk', config);
end

%% plotting error of estimated parameters. Each plot is one parameter which has 5 groups and 3 conditions 
PlotErrorBarPerParam(AllYoungParams, AllHealthyOldParams, AllMCIPosParams, AllMCINegParams, AllMCIUnkParams, config)

%% Two-way ANOVA on (k3-mean return angle) across groups and conditions
%TwowayAnovaOnLastTurn(AllYoungTurnAngleDiff, AllHealthyOldTurnAngleDiff, AllMCIPosTurnAngleDiff, AllMCINegTurnAngleDiff, AllMCIUnkTurnAngleDiff, config);

%% Two-way anova on estmated parameters across 5 groups and 3 conditions
TwowayAnovaOnParams(AllYoungParams, AllHealthyOldParams, AllMCIPosParams, AllMCINegParams, AllMCIUnkParams, config);

%% visulization by merging all MCI
%Merging MCI
AllMCI = MergeMCI(AllMCIPosParams, AllMCINegParams, AllMCIUnkParams);
% performing Two-way anova on MCI v.s. HealthyOld v.s. Young under 3 conditions
TwowayAnovaOnParamsMergeMCI(AllYoungParams, AllHealthyOldParams, AllMCI, config);
% error bar plot of MCI v.s. HealthyOld v.s. Young
PlotErrorBarPerParamMergeMCI(AllYoungParams, AllHealthyOldParams, AllMCI, config)

%% a function for 
function [AllParams, AllX, AllDX, AllTheta, AllflagOoB, AllIC] = getResultsAllConditions(TransformedData, config)
    %get the estimated parameters, X, DX, Theta, IC for all trial
    %conditions for each group of data.
    AllParams = cell(0); AllX = cell(0);AllDX = cell(0); AllTheta = cell(0); AllflagOoB = cell(0); AllIC = cell(0);
    for TRIAL_FILTER=1:3
        config.TrialFilter = TRIAL_FILTER;
        tic
        disp('%%%%%%%%%%%%%%% PERFORMING FITTING %%%%%%%%%%%%%%%');
        Results = PerformGroupFit(TransformedData, config);
%         AllParams{TRIAL_FILTER} = Results.estimatedParams; 
%         AllX{TRIAL_FILTER} = Results.X;
%         AllDX{TRIAL_FILTER}=Results.DX;      
%         AllTheta{TRIAL_FILTER}=Results.THETADX;
%         AllflagOoB{TRIAL_FILTER}=Results.flagOoB;
%         AllIC{TRIAL_FILTER}=Results.IC;
        toc
    end
end    
