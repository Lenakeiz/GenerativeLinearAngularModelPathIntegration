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
savefolder = "C:/Users/Zilong/Desktop/path integration model/Andrea's matlab code/GammaModelAllReplaceWithBias/Output/";

%% setting the configuration
config.UseGlobalSearch = true;
config.USEOoBTrials = true;
%load the data
if config.USEOoBTrials== true 
    YoungControls = TransformPathsOoB(YoungControls); 
else
    YoungControls = TransformPaths(YoungControls);
end

resultfolder = savefolder+"ModelComp/";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% the allocentric PI model without weber's law
config.ModelName = "allocentric";
config.NumParams = 1;
[~, ~, ~, ~, AllYoungIC_Allocentric] = getResultsAllConditions(YoungControls, config);

%% the allocentric PI model with weber's law
config.ModelName = "allocentric_weber";
config.NumParams = 1;
[~, ~, ~, ~, AllYoungIC_AllocentricWeber] = getResultsAllConditions(YoungControls, config);

%% Our base model
%no consider of the rotation gain factor in leg2, i.e., gamma, G3=1, g2=1, g3, k3, sigma, nu. #params=6
config.ModelName = "set_g2_1";
config.NumParams = 5;
[~, ~, ~, ~, AllYoungIC_Base] = getResultsAllConditions(YoungControls, config);

%% Bar Plot of AIC BIC NLL Likelihood
%Plot AIC
ICType = "AIC";
[mean_Allocentric, sem_Allocentric] = getMeanSem(AllYoungIC_Allocentric, ICType);
[mean_AllocentricWeber, sem_AllocentricWeber] = getMeanSem(AllYoungIC_AllocentricWeber, ICType);
[mean_Base, sem_Base] = getMeanSem(AllYoungIC_Base, ICType);
mean_All = [mean_Allocentric',mean_AllocentricWeber', mean_Base'];
sem_All = [sem_Allocentric', sem_AllocentricWeber', sem_Base'];
plotErrorBarFunc(mean_All, sem_All, ICType, resultfolder);

%Plot BIC
ICType = "BIC";
[mean_Allocentric, sem_Allocentric] = getMeanSem(AllYoungIC_Allocentric, ICType);
[mean_AllocentricWeber, sem_AllocentricWeber] = getMeanSem(AllYoungIC_AllocentricWeber, ICType);
[mean_Base, sem_Base] = getMeanSem(AllYoungIC_Base, ICType);
mean_All = [mean_Allocentric',mean_AllocentricWeber',mean_Base'];
sem_All = [sem_Allocentric',sem_AllocentricWeber',sem_Base'];
plotErrorBarFunc(mean_All, sem_All, ICType, resultfolder);

%Plot NEGLL
ICType = "Negll";
[mean_Allocentric, sem_Allocentric] = getMeanSem(AllYoungIC_Allocentric, ICType);
[mean_AllocentricWeber, sem_AllocentricWeber] = getMeanSem(AllYoungIC_AllocentricWeber, ICType);
[mean_Base, sem_Base] = getMeanSem(AllYoungIC_Base, ICType);
mean_All = [mean_Allocentric',mean_AllocentricWeber',mean_Base'];
sem_All = [sem_Allocentric',sem_AllocentricWeber',sem_Base'];
plotErrorBarFunc(mean_All, sem_All, ICType, resultfolder);

%Plot Likelihood
ICType = "Likelihood";
[mean_Allocentric, sem_Allocentric] = getMeanSem(AllYoungIC_Allocentric, ICType);
[mean_AllocentricWeber, sem_AllocentricWeber] = getMeanSem(AllYoungIC_AllocentricWeber, ICType);
[mean_Base, sem_Base] = getMeanSem(AllYoungIC_Base, ICType);
mean_All = [mean_Allocentric',mean_AllocentricWeber',mean_Base'];
sem_All = [sem_Allocentric',sem_AllocentricWeber',sem_Base'];
plotErrorBarFunc(mean_All, sem_All, ICType, resultfolder);

%% A function for getting Results from All Conditions
function [AllParams, AllX, AllDX, AllTheta, AllIC] = getResultsAllConditions(TransformedData, config)
    %get the estimated parameters, X, DX, Theta, IC for all trial
    %conditions for each group of data.
    AllParams = cell(0); AllX = cell(0);AllDX = cell(0); AllTheta = cell(0); AllIC = cell(0);
    for TRIAL_FILTER=1:3
        config.TrialFilter = TRIAL_FILTER;
        tic
        disp('%%%%%%%%%%%%%%% PERFORMING FITTING %%%%%%%%%%%%%%%');
        Results = PerformGroupFit(TransformedData, config);
        AllParams{TRIAL_FILTER} = Results.estimatedParams; 
        AllX{TRIAL_FILTER} = Results.X;
        AllDX{TRIAL_FILTER}=Results.DX;      
        AllTheta{TRIAL_FILTER}=Results.THETADX;
        AllIC{TRIAL_FILTER}=Results.IC;
        toc
    end
end    

function [mean_ic, sem_ic] = getMeanSem(IC, ICType)
    mean_ic = zeros(1,3); sem_ic = zeros(1,3);
    for i=1:3 %three conditions
        if ICType == "AIC"
            aicdata = [];
            for j=1:length(IC{i})
                aicdata = [aicdata,IC{i}{j}.aic];
            end
            mean_ic(i) = mean(aicdata);
            sem_ic(i) = std(aicdata)./sqrt(length(aicdata));

        elseif ICType == "BIC"
            bicdata = [];
            for j=1:length(IC{i})
                bicdata = [bicdata,IC{i}{j}.bic];
            end
            mean_ic(i) = mean(bicdata);
            sem_ic(i) = std(bicdata)./sqrt(length(bicdata));

        elseif ICType == "Negll"
            neglldata = [];
            for j=1:length(IC{i})
                neglldata = [neglldata,IC{i}{j}.negll];
            end
            mean_ic(i) = mean(neglldata);
            sem_ic(i) = std(neglldata)./sqrt(length(neglldata));
        elseif ICType == "Likelihood"
            likelihooddata = [];
            for j=1:length(IC{i})
                likelihooddata = [likelihooddata,IC{i}{j}.likelihood];
            end
            mean_ic(i) = mean(likelihooddata);
            sem_ic(i) = std(likelihooddata)./sqrt(length(likelihooddata));            
        else
            error("Choose correct IC type!")
        end
    end
end

function plotErrorBarFunc(mean_all, sem_all, ICType, resultfolder)
    
    f = figure('visible','off','Position', [100 100 500 300]);

    b = bar(mean_all, 'grouped');
    hold on
    
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(mean_all);
    % Get the x coordinate of the bars
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end
    % Plot the errorbars
    errorbar(x',mean_all,sem_all,'k','linestyle','none');
    
    legend(b, {'allocentric', 'allocentric-weber', 'base'}, 'Location','northwest','FontSize',12, 'NumColumns',3);
    
    if ICType=="AIC"
        ylabel('Akaike information criterion (AIC)', 'FontName', 'Arial', 'FontSize',12);
    elseif ICType=="BIC"
        ylabel('Bayesian Inference Criterion (BIC)', 'FontName', 'Arial', 'FontSize',12);
    elseif ICType=="Negll"
        ylabel('Negative Loglikelihood','FontName', 'Arial', 'FontSize',12);
    elseif ICType=="Likelihood"
        ylabel('Likelihood','FontName', 'Arial', 'FontSize',12);
    end
    
    %ylim([4,7]);

    %xticklabels({'No change', 'No distal cue', 'No optical flow'});
    set(gca,'XTickLabel',{'No change' 'No distal cue' 'No optical flow'}, 'FontName','Arial','FontSize',12);
    exportgraphics(f,resultfolder+ICType+".png",'Resolution',300);
end
