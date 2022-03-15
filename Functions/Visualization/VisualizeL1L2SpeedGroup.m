%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default');

disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');

%%
config.Speed.alpha = 0.9; %Paramanter for running speed calculation
config.Speed.timeOffsetAfterFlagReach = 1.5; %Time to track after flag reached in seconds 
config.Speed.smoothWindow = 10; %a trarate should be 10Hz so 4 secs window is 40 datapoints
config.Speed.velocityCutoff = 0.2; % velocity cutoff to select only the walking part of the reconstructed velocity
config.Speed.timeOffsetForDetectedTemporalWindow = 0.4; % time in seconds that will push earlier/ the detected rising edge
config.Speed.tresholdForBadParticipantL1Recontruction = 2.0; % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
%%
config.Speed.tresholdForBadParticipantL1Recontruction = 1.55;
disp('%%%%%%%%%%%%%%% Calculating Tracking Path - YOUNG ... %%%%%%%%%%%%%%%');
YoungControls   = CalculateTrackingPath(YoungControls, config);
config.Speed.tresholdForBadParticipantL1Recontruction = 2.0;
disp('%%%%%%%%%%%%%%% Calculating Tracking Path - HEALTHYCONTROLS ... %%%%%%%%%%%%%%%');
HealthyControls = CalculateTrackingPath(HealthyControls, config);
disp('%%%%%%%%%%%%%%% Calculating Tracking Path - UNKNOWN ... %%%%%%%%%%%%%%%');
Unknown         = CalculateTrackingPath(Unknown, config);
disp('%%%%%%%%%%%%%%% Calculating Tracking Path - MCI POS ... %%%%%%%%%%%%%%%');
MCIPos          = CalculateTrackingPath(MCIPos, config);
disp('%%%%%%%%%%%%%%% Calculating Tracking Path - MCI NEG ... %%%%%%%%%%%%%%%');
MCINeg          = CalculateTrackingPath(MCINeg, config);
ColorPattern;
disp('%%%%%%%%%%%%%%% DONE CALCULATING TRACKING PATH ... %%%%%%%%%%%%%%%');
%%
disp('%%%%%%%%%%%%%%% Plotting smoothed velocity %%%%%%%%%%%%%%%');
plotCumulativeSpeedL1L2(YoungControls,'YoungControls',config);
plotCumulativeSpeedL1L2(HealthyControls,'HealthyControls',config);
plotCumulativeSpeedL1L2(Unknown,'Unknown',config);
plotCumulativeSpeedL1L2(MCINeg,'MCINeg',config);
plotCumulativeSpeedL1L2(MCIPos,'MCIPos',config);
disp('%%%%%%%%%%%%%%% Plotting smoothed velocity - FINISHED %%%%%%%%%%%%%%%');
%%
disp('%%%%%%%%%%%%%%% Plotting reconstructed lenghts SMOOTHED %%%%%%%%%%%%%%%');
plotGroupwiseLenghtReconstructionForL1L2Smoothed(YoungControls,"Young controls",config);
plotGroupwiseLenghtReconstructionForL1L2Smoothed(HealthyControls,"Healthy controls",config);
plotGroupwiseLenghtReconstructionForL1L2Smoothed(Unknown,"MCI Unknown",config);
plotGroupwiseLenghtReconstructionForL1L2Smoothed(MCINeg,"MCI Negative",config);
plotGroupwiseLenghtReconstructionForL1L2Smoothed(MCIPos,"MCI Positive",config);
disp('%%%%%%%%%%%%%%% Plotting reconstructed lenghts SMOOTHED - FINISHED %%%%%%%%%%%%%%%');

%%
disp('%%%%%%%%%%%%%%% Plotting reconstructed lenghts UNSMOOTHED %%%%%%%%%%%%%%%');
plotGroupwiseLenghtReconstructionForL1L2Unsmoothed(YoungControls,"Young controls",config);
plotGroupwiseLenghtReconstructionForL1L2Unsmoothed(HealthyControls,"Healthy controls",config);
plotGroupwiseLenghtReconstructionForL1L2Unsmoothed(Unknown,"MCI Unknown",config);
plotGroupwiseLenghtReconstructionForL1L2Unsmoothed(MCINeg,"MCI Negative",config);
plotGroupwiseLenghtReconstructionForL1L2Unsmoothed(MCIPos,"MCI Positive",config);
disp('%%%%%%%%%%%%%%% Plotting reconstructed lenghts UNSMOOTHED - FINISHED %%%%%%%%%%%%%%%');

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
function plotGroupwiseLenghtReconstructionForL1L2Smoothed(Group,GroupName,config)
    
    savefolder = pwd + "/Output/";
    % setting the configuration
    resultfolder = savefolder+"PaperFigs/ReconstructedL1L2Smoothed/";
    config.ResultFolder = resultfolder;
    %create storing folder for trajectory if not exist
    if ~exist(resultfolder, 'dir')
       mkdir(resultfolder);
    end

    allData = vertcat(Group.Reconstructed{:});
    allData = rmmissing(allData);

    jitter_value = 0.2*(rand(height(allData),1)-0.5);

    % set figure info
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
    pt = plot([ [ones(height(allData),1)+jitter_value]'; [2*ones(height(allData),1) + jitter_value]'; [3*ones(height(allData),1) + jitter_value]' ],...
                [allData.L1Real'; allData.L1Smoothed'; allData.L1SmoothedFiltered'],...
                '-', 'Color', [config.color_scheme_npg(2,:) 0.2],'linewidth',2);

    st = scatter( [ones(height(allData),1)+jitter_value, 2*ones(height(allData),1) + jitter_value, 3*ones(height(allData),1) + jitter_value],...
                  [allData.L1Real, allData.L1Smoothed, allData.L1SmoothedFiltered], 50,...
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
    pt = plot([ [ones(height(allData),1)+jitter_value]'; [2*ones(height(allData),1) + jitter_value]'; [3*ones(height(allData),1) + jitter_value]' ],...
                [allData.L2Real'; allData.L2Smoothed'; allData.L2SmoothedFiltered'],...
                '-', 'Color', [config.color_scheme_npg(2,:) 0.2],'linewidth',2);

    st = scatter( [ones(height(allData),1)+jitter_value, 2*ones(height(allData),1) + jitter_value, 3*ones(height(allData),1) + jitter_value],...
                  [allData.L2Real, allData.L2Smoothed, allData.L2SmoothedFiltered], 50,...
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

%% Plotting recontructed lenght over speed with unfiltered data (bear in mind however that the parameters for this are extracted from the filtered data)
function plotGroupwiseLenghtReconstructionForL1L2Unsmoothed(Group,GroupName,config)
    
    savefolder = pwd + "/Output/";
    % setting the configuration
    resultfolder = savefolder+"PaperFigs/ReconstructedL1L2Unsmoothed/";
    config.ResultFolder = resultfolder;
    %create storing folder for trajectory if not exist
    if ~exist(resultfolder, 'dir')
       mkdir(resultfolder);
    end

    allData = vertcat(Group.Reconstructed{:});
    allData = rmmissing(allData);

    jitter_value = 0.2*(rand(height(allData),1)-0.5);

    % set figure info
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
    pt = plot([ [ones(height(allData),1)+jitter_value]'; [2*ones(height(allData),1) + jitter_value]'; [3*ones(height(allData),1) + jitter_value]' ],...
                [allData.L1Real'; allData.L1Unsmoothed'; allData.L1UnsmoothedFiltered'],...
                '-', 'Color', [config.color_scheme_npg(2,:) 0.2],'linewidth',2);

    st = scatter( [ones(height(allData),1)+jitter_value, 2*ones(height(allData),1) + jitter_value, 3*ones(height(allData),1) + jitter_value],...
                  [allData.L1Real, allData.L1Unsmoothed, allData.L1UnsmoothedFiltered], 50,...
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
    pt = plot([ [ones(height(allData),1)+jitter_value]'; [2*ones(height(allData),1) + jitter_value]'; [3*ones(height(allData),1) + jitter_value]' ],...
                [allData.L2Real'; allData.L2Unsmoothed'; allData.L2UnsmoothedFiltered'],...
                '-', 'Color', [config.color_scheme_npg(2,:) 0.2],'linewidth',2);

    st = scatter( [ones(height(allData),1)+jitter_value, 2*ones(height(allData),1) + jitter_value, 3*ones(height(allData),1) + jitter_value],...
                  [allData.L2Real, allData.L2Unsmoothed, allData.L2UnsmoothedFiltered], 50,...
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

%% Plotting the distribution of times for each of the participants