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
%Model_Name = "base"; numParams=7;   

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

config.USEOoBTrials = true;
if config.USEOoBTrials==true
    resultfolder = savefolder+"ModelComparison/OoB/";
else
    resultfolder = savefolder+Model_Name+"ModelComparison/NoOoB/";
end
config.ResultFolder = resultfolder;

config.TrialFilter = 3; %only use the no change condition trials
%create storing folder for store data if not exist
config.TrialFolder = resultfolder+"TrailType"+num2str(config.TrialFilter)+"/";
if ~exist(config.TrialFolder, 'dir')
   mkdir(config.TrialFolder);
end

%plot
config.ifReturnAngleDistribution = true;
config.ifPhysicalTrajectory = true;
config.ifMentalTrajectory = true;
config.ifEndpointDistribution = true;


%% Transform the data
if config.USEOoBTrials==true
    YoungControls = TransformPathsOoB(YoungControls);
else
    YoungControls = TransformPaths(YoungControls);
end


%% Performing data fitting on the base_add_G3 model
[AIC1, BIC1, NEGLL1] = getICOfModel("base_add_G3", 8, YoungControls, config);

%% Performing data fitting on the base model
[AIC2, BIC2, NEGLL2] = getICOfModel("base", 7, YoungControls, config);

%% Performing data fitting on the set_g3_1_k3_0 model
[AIC3, BIC3, NEGLL3] = getICOfModel("set_g3_1_k3_0", 5, YoungControls, config);

%% Performing data fitting on the set_g2_1 model
[AIC4, BIC4, NEGLL4] = getICOfModel("set_g2_1", 6, YoungControls, config);

%% Performing data fitting on the set_g2_1_g3_1_k3_0 model
[AIC5, BIC5, NEGLL5] = getICOfModel("set_g2_1_g3_1_k3_0", 4, YoungControls, config);

%% Performing data fitting on the set_bG2_1 model
[AIC6, BIC6, NEGLL6] = getICOfModel("set_bG2_1", 6, YoungControls, config);

%% Performing data fitting on the set_bG1_1 model
[AIC7, BIC7, NEGLL7] = getICOfModel("set_bG1_1", 6, YoungControls, config);

%% Performing data fitting on the set_same_bG model
[AIC8, BIC8, NEGLL8] = getICOfModel("set_same_bG", 6, YoungControls, config);

%% Performing data fitting on the set_same_sg model
[AIC9, BIC9, NEGLL9] = getICOfModel("set_same_sg", 5, YoungControls, config);

%% Performing data fitting on the set_same_Gg model
[AIC10, BIC10, NEGLL10] = getICOfModel("set_same_Gg", 4, YoungControls, config);

%% plot AIC
mean_all = [AIC1.mean,AIC2.mean,AIC3.mean,AIC4.mean,AIC5.mean,AIC6.mean,AIC7.mean,AIC8.mean,AIC9.mean,AIC10.mean];
sem_all = [AIC1.sem,AIC2.sem,AIC3.sem,AIC4.sem,AIC5.sem,AIC6.sem,AIC7.sem,AIC8.sem,AIC9.sem,AIC10.sem];
plotErrorBarFunc(mean_all, sem_all, "aic", config)

%% plot BIC
mean_all = [BIC1.mean,BIC2.mean,BIC3.mean,BIC4.mean,BIC5.mean,BIC6.mean,BIC7.mean,BIC8.mean,BIC9.mean,BIC10.mean];
sem_all = [BIC1.sem,BIC2.sem,BIC3.sem,BIC4.sem,BIC5.sem,BIC6.sem,BIC7.sem,BIC8.sem,BIC9.sem,BIC10.sem];
plotErrorBarFunc(mean_all, sem_all, "bic", config)

%% plot NEGLL
mean_all = [NEGLL1.mean,NEGLL2.mean,NEGLL3.mean,NEGLL4.mean,NEGLL5.mean,NEGLL6.mean,NEGLL7.mean,NEGLL8.mean,NEGLL9.mean,NEGLL10.mean];
sem_all = [NEGLL1.sem,NEGLL2.sem,NEGLL3.sem,NEGLL4.sem,NEGLL5.sem,NEGLL6.sem,NEGLL7.sem,NEGLL8.sem,NEGLL9.sem,NEGLL10.sem];
plotErrorBarFunc(mean_all, sem_all, "negll", config)

%%
function plotErrorBarFunc(mean_all, sem_all, ICType, config)

    f = figure('visible','off','Position', [100 100 1000 800]);
    
    b = bar(mean_all);
    hold on
    
    % Calculate the number of groups and number of bars in each group
    nbars = length(mean_all);
    % Get the x coordinate of the bars
    x = 1:nbars;
    % Plot the errorbars
    errorbar(x,mean_all,sem_all,'k','linestyle','none');
    
    %ylabel('Bayesian Inference Criterion (BIC)', 'FontSize',20);
    %ylabel('Likelihood', 'FontSize',20);
    ylabel(ICType, 'FontSize',20);
    title(ICType, 'FontSize',20);
    
    %ylim([2,4]);
    
    set(gca, 'XTickLabel', {'Base-add-G3' ...
                            'Base' 'g_3=1 k_3=0' ...
                            'g_2=1' 'g_2=1 g_3=1 k_3=0' ...
                            'G_2=1' 'G_1=1'  'Same-bg' ...
                            'Same-sg' 'Same-Gg'});
    
    exportgraphics(f,config.TrialFolder+ICType+".png",'Resolution',300);

end
   
%%
function [AIC, BIC, NEGLL] = getICOfModel(ModelName, NumParams, YoungControls, config)
    config.ModelName = ModelName;
    config.NumParams = NumParams;
    Results = PerformGroupFit(YoungControls, config);
    IC = Results.IC;

    AICInfo = zeros(1,length(IC)); BICInfo=zeros(1,length(IC)); NEGLLInfo=zeros(1,length(IC));
    for i=1:length(IC)
        AICInfo(i) = IC{i}.aic;
        BICInfo(i) = IC{i}.bic;
        NEGLLInfo(i) = IC{i}.negll;
    end

    AIC.mean = mean(AICInfo); AIC.sem = std(AICInfo)/sqrt(length(AICInfo));
    BIC.mean = mean(BICInfo); BIC.sem = std(BICInfo)/sqrt(length(BICInfo));
    NEGLL.mean = mean(NEGLLInfo); NEGLL.sem = std(NEGLLInfo)/sqrt(length(NEGLLInfo));
end