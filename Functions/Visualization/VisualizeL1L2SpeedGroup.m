%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default');

disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');

%%
config.Speed.alpha = 0.9; %Paramanter for running speed calculation
config.Speed.timeOffsetAfterFlagReach = 0; %Time to track after flag reached in seconds 
config.Speed.smoothWindow = 10; % trarate should be 10Hz so 4 secs window is 40 datapoints
config.Speed.velocityCutoff = 0.2; % velocity cutoff to select only the walking part of the reconstructed velocity
config.Speed.timeOffsetForDetectedTemporalWindow = 0.5; % time in seconds that will push earlier/ the detected rising edge
%%
disp('%%%%%%%%%%%%%%% Calculating Tracking Path - YOUNG ... %%%%%%%%%%%%%%%');
YoungControls   = CalculateTrackingPath(YoungControls, config);
disp('%%%%%%%%%%%%%%% Calculating Tracking Path - HEALTHYCONTROLS ... %%%%%%%%%%%%%%%');
HealthyControls = CalculateTrackingPath(HealthyControls, config);
disp('%%%%%%%%%%%%%%% Calculating Tracking Path - UNKNOWN ... %%%%%%%%%%%%%%%');
Unknown         = CalculateTrackingPath(Unknown, config);
disp('%%%%%%%%%%%%%%% Calculating Tracking Path - MCI POS ... %%%%%%%%%%%%%%%');
MCIPos          = CalculateTrackingPath(MCIPos, config);
disp('%%%%%%%%%%%%%%% Calculating Tracking Path - MCI NEG ... %%%%%%%%%%%%%%%');
MCINeg          = CalculateTrackingPath(MCINeg, config);
disp('%%%%%%%%%%%%%%% DONE CALCULATING TRACKING PATH ... %%%%%%%%%%%%%%%');
%%
ColorPattern; 
%%
disp('%%%%%%%%%%%%%%% Plotting smoothed velocity %%%%%%%%%%%%%%%');
plotCumulativeSpeedL1L2(YoungControls,'YoungControls',config);
plotCumulativeSpeedL1L2(HealthyControls,'HealthyControls',config);
plotCumulativeSpeedL1L2(Unknown,'Unknown',config);
plotCumulativeSpeedL1L2(MCINeg,'MCINeg',config);
%%
plotCumulativeSpeedL1L2(MCIPos,'MCIPos',config);
disp('%%%%%%%%%%%%%%% Plotting smoothed velocity - FINISHED %%%%%%%%%%%%%%%');

%%
disp('%%%%%%%%%%%%%%% Plotting reconstructed lenghts %%%%%%%%%%%%%%%');
youngBadPpt           = plotGroupwiseLenghtReconstructionForL1L2(YoungControls,"Young controls",config,1.55);
healthyControlsBadPpt = plotGroupwiseLenghtReconstructionForL1L2(HealthyControls,"Healthy controls",config,2.0);
unknownBadPpt         = plotGroupwiseLenghtReconstructionForL1L2(Unknown,"MCI Unknown",config,2.0);
mciNegBadPpt          = plotGroupwiseLenghtReconstructionForL1L2(MCINeg,"MCI Negative",config,2.0);
mciPosBadPpt          = plotGroupwiseLenghtReconstructionForL1L2(MCIPos,"MCI Positive",config,2.0);
disp('%%%%%%%%%%%%%%% Plotting reconstructed lenghts - FINISHED %%%%%%%%%%%%%%%');

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
        
        f = figure('visible','off','Position', [0 0 1200 700]);
        
        %%% Font type and size setting %%%
        % Using Arial as default because all journals normally require the font to
        % be either Arial or Helvetica
        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)
        
        for tr = 1:nTrials
        
            subplot("Position", [0.1 0.1 0.35 0.7]);
            hold on;
        
            plot(Group.TrackedL1{1,pId}{tr,1}.ShiftedTime, Group.TrackedL1{1,pId}{tr,1}.Smoothed_Vel_proj,...Color=[config.color_scheme_npg(9,:) 1],...
                LineStyle="-",...
                LineWidth=1.3);
        
            xlabel('Time (s)');
            ylabel('$\bar{V} \cdot \hat{L_1}$','Interpreter','latex');
        
            title("L_1",FontSize=14);

            ylim([-0.2 1.0]);

            subplot("Position", [0.55 0.1 0.35 0.7]);
            hold on;
        
            plot(Group.TrackedL2{1,pId}{tr,1}.ShiftedTime, Group.TrackedL2{1,pId}{tr,1}.Smoothed_Vel_proj,...Color=[config.color_scheme_npg(9,:) 1],...
                LineStyle="-",...
                LineWidth=1.3);
        
            xlabel('Time (s)');
            ylabel('$\bar{V} \cdot \hat{L_2}$','Interpreter','latex');
        
            title("L_2",FontSize=14);

            ylim([-0.2 1.0]);

        end
        
        hold off
        
        sTitle = sgtitle("Gaussian filtered reconstructed speed - pID " + num2str(pId));
        sTitle.FontSize = 15;
        sTitle.FontWeight = 'bold';
        
        exportgraphics(f,config.ResultFolder+"/pID_"+num2str(pId)+"_SpeedL1L2"+".png",'Resolution',300);

    end

end

%% Plotting recontructed lenght over speed
function [groupFaultyIndices] = plotGroupwiseLenghtReconstructionForL1L2(Group,GroupName,config,thresholdForFaultyParticipants)
    
    savefolder = pwd + "/Output/";
    % setting the configuration
    resultfolder = savefolder+"PaperFigs/ReconstructedL1L2/";
    config.ResultFolder = resultfolder;
    %create storing folder for trajectory if not exist
    if ~exist(resultfolder, 'dir')
       mkdir(resultfolder);
    end

    groupFaultyIndices = [];

    pSize = size(Group.TrackedL1,2);
    l1reconstructed = {};
    l1reconstructedFilter = {};
    l1real = {};

    l2reconstructed = {};
    l2reconstructedFilter = {};
    l2real = {};
    
    for pId = 1:pSize

        trialSize = size(Group.TrackedL1{1,pId},1);

        for trialId = 1:trialSize

            Cone_pos    = Group.FlagPos{1,pId}{trialId,1};
            currl1real = norm(Cone_pos(2,[1 3]) - Cone_pos(1,[1 3]));
            currl2real = norm(Cone_pos(3,[1 3]) - Cone_pos(2,[1 3]));
            
            TrackedL1Smoothed = Group.TrackedL1{1,pId}{trialId,1}.Smoothed_Vel_proj;
            currL1rec         = trapz(Group.TrackedL1{1,pId}{trialId,1}.Time,TrackedL1Smoothed);
            currL1recfilter   = trapz(Group.TrackedL1{1,pId}{trialId,1}.Time(Group.TrackedL1{1,pId}{trialId,1}.Filthered_Vel_proj),TrackedL1Smoothed(Group.TrackedL1{1,pId}{trialId,1}.Filthered_Vel_proj));

            TrackedL2Smoothed = Group.TrackedL2{1,pId}{trialId,1}.Smoothed_Vel_proj;
            currL2rec         = trapz(Group.TrackedL2{1,pId}{trialId,1}.Time,TrackedL2Smoothed);
            currL2recfilter   = trapz(Group.TrackedL2{1,pId}{trialId,1}.Time(Group.TrackedL2{1,pId}{trialId,1}.Filthered_Vel_proj),TrackedL2Smoothed(Group.TrackedL2{1,pId}{trialId,1}.Filthered_Vel_proj));

            if(currl1real < thresholdForFaultyParticipants)
                groupFaultyIndices = [groupFaultyIndices, pId];
                l1real{pId,1}{trialId,1} = nan;
                l1reconstructed{pId,1}{trialId,1} = nan;
                l1reconstructedFilter{pId,1}{trialId,1} = nan;
    
                l2real{pId,1}{trialId,1} = nan;
                l2reconstructed{pId,1}{trialId,1} = nan;
                l2reconstructedFilter{pId,1}{trialId,1} = nan;
            else
                l1real{pId,1}{trialId,1} = currl1real;
                l1reconstructed{pId,1}{trialId,1} = currL1rec;
                l1reconstructedFilter{pId,1}{trialId,1} = currL1recfilter;
                
                l2real{pId,1}{trialId,1} = currl2real;
                l2reconstructed{pId,1}{trialId,1} = currL2rec;
                l2reconstructedFilter{pId,1}{trialId,1} = currL2recfilter;
            end

        end
    
    end

    groupFaultyIndices = unique(groupFaultyIndices);

    l1real          = cell2mat(vertcat(l1real{:}));
    l1reconstructed = cell2mat(vertcat(l1reconstructed{:}));
    l1reconstructedFilter = cell2mat(vertcat(l1reconstructedFilter{:}));


    l2reconstructed = cell2mat(vertcat(l2reconstructed{:}));
    l2real          = cell2mat(vertcat(l2real{:}));    
    l2reconstructedFilter = cell2mat(vertcat(l2reconstructedFilter{:}));

    jitter_value = 0.2*(rand(height(l1reconstructed),1)-0.5);

    %% set figure info
    f = figure('visible','off','Position', [0 0 1200 700]);
    %f = figure('Position', [0 0 1400 600]);

    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)

    %%% Color definition %%%
    %pairline_color = 'k';

    ax_left = subplot("Position", [0.1 0.1 0.35 0.7]);

    hold on

    % Tutte le x dei punti a sinistra come prima riga, tutte le x dei punti a destra come seconda riga 
    pt = plot([ [ones(height(l1reconstructed),1)+jitter_value]'; [2*ones(height(l1reconstructed),1) + jitter_value]'; [3*ones(height(l1reconstructed),1) + jitter_value]' ],...
                [l1real'; l1reconstructed'; l1reconstructedFilter'],...
                '-', 'Color', [config.color_scheme_npg(2,:) 0.2],'linewidth',2);

    st = scatter( [ones(height(l1reconstructed),1)+jitter_value, 2*ones(height(l1reconstructed),1) + jitter_value, 3*ones(height(l1reconstructed),1) + jitter_value],...
                  [l1real, l1reconstructed, l1reconstructedFilter], 50,...
                  "filled","MarkerEdgeColor","k","MarkerFaceColor",config.color_scheme_npg(1,:),...
                  "MarkerEdgeAlpha",0.8,"MarkerFaceAlpha",0.3);

    hold off    

    ax_left.XAxis.TickValues = [1 2 3];
    ax_left.XAxis.TickLabels = {'Real', 'Reconstructed', 'Filtered'};

    title("L_{1}",FontSize=14)

    ylim([1.0 6.0]);

    ax_right = subplot("Position", [0.55 0.1 0.35 0.7]);
    ax_right.FontSize = 12;

    hold on

    % Tutte le x dei punti a sinistra come prima riga, tutte le x dei punti a destra come seconda riga 
    pt = plot([ [ones(height(l2reconstructed),1)+jitter_value]'; [2*ones(height(l2reconstructed),1) + jitter_value]'; [3*ones(height(l2reconstructed),1) + jitter_value]' ],...
                [l2real'; l2reconstructed'; l2reconstructedFilter'],...
                '-', 'Color', [config.color_scheme_npg(2,:) 0.2],'linewidth',2);

    st = scatter( [ones(height(l2reconstructed),1)+jitter_value, 2*ones(height(l2reconstructed),1) + jitter_value, 3*ones(height(l2reconstructed),1) + jitter_value],...
                  [l2real, l2reconstructed, l2reconstructedFilter], 50,...
                  "filled","MarkerEdgeColor","k","MarkerFaceColor",config.color_scheme_npg(1,:),...
                  "MarkerEdgeAlpha",0.8,"MarkerFaceAlpha",0.3);

    hold off

    title("L_2",FontSize=14)

    ax_right.XAxis.TickValues = [1 2 3];
    ax_right.XAxis.TickLabels = {'Real', 'Reconstructed', 'Filtered'};

    ylim([1.0 6.0]);

    sTitle = sgtitle("Real vs Reconstructed Outbound Path - " + GroupName);
    sTitle.FontSize = 15;
    sTitle.FontWeight = 'bold';

    exportgraphics(f,config.ResultFolder+"/L1L2_reconstruction_"+ GroupName + ".png",'Resolution',300);

end