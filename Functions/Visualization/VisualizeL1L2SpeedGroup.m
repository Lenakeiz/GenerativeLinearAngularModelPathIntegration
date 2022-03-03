%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default');

disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');

%%
config.Speed.alpha = 0.9; %Paramanter for running speed calculation
config.Speed.TOffsetAfterFlagReach = 2; %Time to track after flag reached in seconds 
config.Speed.smoothWindow = 40; % tracking rate should be 10Hz so 4 secs window

%%
YoungControls   = CalculateTrackingPath(YoungControls, config);
HealthyControls = CalculateTrackingPath(HealthyControls, config);
Unknown         = CalculateTrackingPath(Unknown, config);
MCIPos          = CalculateTrackingPath(MCIPos, config);
MCINeg          = CalculateTrackingPath(MCINeg, config);

plotCumulativeSpeedL1L2(YoungControls,'YoungControls',config);
plotCumulativeSpeedL1L2(HealthyControls,'HealthyControls',config);
plotCumulativeSpeedL1L2(Unknown,'Unknown',config);
plotCumulativeSpeedL1L2(MCIPos,'MCIPos',config);
plotCumulativeSpeedL1L2(MCINeg,'MCINeg',config);

%% Plotting cumulative plots
function plotCumulativeSpeedL1L2(Group,GroupName,config)
    
    savefolder = pwd + "/Output/";
    % setting the configuration
    resultfolder = savefolder+"PaperFigs/GaussianFilteredSpeedL1L2/" + GroupName;
    config.ResultFolder = resultfolder;
    %create storing folder for trajectory if not exist
    if ~exist(resultfolder, 'dir')
       mkdir(resultfolder);
    end
   
    nParticipants = size(Group.TrackedL1,2);
    
    for pId = 1:nParticipants

        nTrials = size(Group.TrackedL1{1,pId},1);
        
        f = figure('visible','off','Position', [0 0 1400 600]);
        
        %%% Font type and size setting %%%
        % Using Arial as default because all journals normally require the font to
        % be either Arial or Helvetica
        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)
        
        for tr = 1:nTrials
        
            subplot("Position", [0.05 0.05 0.4 0.8]);
            hold on;
            maxL1 = max(Group.TrackedL1{1,pId}{tr,1}.Smoothed_Vel_proj);
            idxMaxL1 = find(Group.TrackedL1{1,pId}{tr,1}.Smoothed_Vel_proj == maxL1);
            
            tempTime = Group.TrackedL1{1,pId}{tr,1}.Time - Group.TrackedL1{1,pId}{tr,1}.Time(idxMaxL1);
        
            plot(tempTime, Group.TrackedL1{1,pId}{tr,1}.Smoothed_Vel_proj,...Color=[config.color_scheme_npg(9,:) 1],...
                LineStyle="-",...
                LineWidth=1.3);
        
            xlabel('Time (s)');
            ylabel('$\bar{V} \cdot \hat{L1}$','Interpreter','latex');
        
            subplot("Position", [0.55 0.05 0.4 0.8]);
            hold on;
        
            maxL2 = max(Group.TrackedL2{1,pId}{tr,1}.Smoothed_Vel_proj);
            idxMaxL2 = find(Group.TrackedL2{1,pId}{tr,1}.Smoothed_Vel_proj == maxL2);
            
            tempTime = Group.TrackedL2{1,pId}{tr,1}.Time - Group.TrackedL2{1,pId}{tr,1}.Time(idxMaxL2);
        
            plot(tempTime, Group.TrackedL2{1,pId}{tr,1}.Smoothed_Vel_proj,...Color=[config.color_scheme_npg(9,:) 1],...
                LineStyle="-",...
                LineWidth=1.3);
        
            xlabel('Time (s)');
            ylabel('$\bar{V} \cdot \hat{L2}$','Interpreter','latex');
         
        end
        
        hold off
        
        sTitle = sgtitle("Gaussian filtered reconstructed speed - pID " + num2str(pId));
        sTitle.FontSize = 15;
        sTitle.FontWeight = 'bold';
        
        exportgraphics(f,config.ResultFolder+"/pID_"+num2str(pId)+"_SpeedL1L2"+".png",'Resolution',300);

    end


end