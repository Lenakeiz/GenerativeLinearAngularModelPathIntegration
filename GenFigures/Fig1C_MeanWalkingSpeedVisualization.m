%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default'); %for code reproducibility

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
savefolder = pwd + "/Output/";

%% setting the configuration
config.Speed.alpha                                  = 0.9;                 %Paramanter for running speed calculation
config.Speed.timeOffsetAfterFlagReach               = 1.5;                 %Time to track after flag reached in seconds 
config.Speed.smoothWindow                           = 10;                  % tracking rate should be 10Hz so 4 secs window is 40 datapoints
config.Speed.velocityCutoff                         = 0.2;                 % velocity cutoff to select only the walking part of the reconstructed velocity
config.Speed.timeOffsetForDetectedTemporalWindow    = 0.4;                 % time in seconds that will push earlier/ the detected rising edge
config.UseGlobalSearch = true;
resultfolder = savefolder+"PaperFigs/Fig1C";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% Model fitting for YoungControl data
%% calculating tracking path and transoform data
config.Speed.tresholdForBadParticipantL1Recontruction = 1.55;   % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
YoungControls   = CalculateTrackingPath(YoungControls, config);
YoungControls   = TransformPaths(YoungControls);%transform data
YoungmeanS      = getMeanWalkingSpeed(YoungControls);

%% Model fitting for HealthyOld data
config.Speed.tresholdForBadParticipantL1Recontruction = 2.0; 
HealthyControls = CalculateTrackingPath(HealthyControls, config);
HealthyOld      = TransformPaths(HealthyControls);%transform data
HealthyOldmeanS = getMeanWalkingSpeed(HealthyOld);

%% Model fitting for MCIPos
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCIPos          = CalculateTrackingPath(MCIPos, config);
MCIPos          = TransformPaths(MCIPos);%transform data
MCIPosmeanS     = getMeanWalkingSpeed(MCIPos);

%% Model fitting for MCINeg
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCINeg          = CalculateTrackingPath(MCINeg, config);
MCINeg          = TransformPaths(MCINeg);%transform data
MCINegmeanS     = getMeanWalkingSpeed(MCINeg);

%% Model fitting for MCIUnk
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
Unknown         = CalculateTrackingPath(Unknown, config);
MCIUnk          = TransformPaths(Unknown);%transform data
MCIUnkmeanS     = getMeanWalkingSpeed(MCIUnk);

%% Setting colors for using in plots
ColorPattern; 

%% Two-way anova on mean speed of 5 groups and 3 conditions
config.type = "meanspeed";
[dist_anova_tab,~,~, ~] = TwowayAnovaOnRealData(YoungmeanS, HealthyOldmeanS, MCIPosmeanS, MCINegmeanS, MCIUnkmeanS, config);

%% aBarScatter plot of return distance of all groups and conditions
plotBarScatterOfMeanSpeed(YoungmeanS, HealthyOldmeanS, MCIPosmeanS, MCINegmeanS, MCIUnkmeanS, dist_anova_tab, config)

%% getting Mean Walking Speed of all participants in a group
function meanS = getMeanWalkingSpeed(Dat)
    numSubjs = length(Dat.Reconstructed);
    meanS = NaN(3,numSubjs);
    for subj=1:numSubjs
        
        %filter out participants who did short walking
        if ismember(subj, Dat.BadPptIdxs)
            continue
        end
        
        %condition index
        IdxCond1    = Dat.CondTable{subj}.Condition==1;
        IdxCond2    = Dat.CondTable{subj}.Condition==2;
        IdxCond3    = Dat.CondTable{subj}.Condition==3;

        %length on leg1
        length1     = table2array(Dat.Reconstructed{subj}(:,'L1Real'));
        %duration on leg1
        duration1   = table2array(Dat.Reconstructed{subj}(:,'T_L1'));
        meanspeed1  = length1./duration1;

        %length on leg2
        length2     = table2array(Dat.Reconstructed{subj}(:,'L2Real'));
        %duration on leg2
        duration2   = table2array(Dat.Reconstructed{subj}(:,'T_L2'));
        meanspeed2  = length2./duration2;

        %mean speed of condition 1,2,3
        meanspeed_cond1 = mean([meanspeed1(IdxCond1);meanspeed2(IdxCond1)],'omitnan');
        meanspeed_cond2 = mean([meanspeed1(IdxCond2);meanspeed2(IdxCond2)],'omitnan');
        meanspeed_cond3 = mean([meanspeed1(IdxCond3);meanspeed2(IdxCond3)],'omitnan');
        
        meanS(1,subj) = meanspeed_cond1;
        meanS(2,subj) = meanspeed_cond2;
        meanS(3,subj) = meanspeed_cond3;
    end

    %delete NANs in meanS
    nanacol = isnan(meanS(1,:));
    meanS   = meanS(:,~nanacol);
end

%% plot Bar scatter of Return distance
function plotBarScatterOfMeanSpeed(YoungmeanS, HealthyOldmeanS, MCIPosmeanS, MCINegmeanS, MCIUnkmeanS, dist_anova_tab, config)
    % YoungmeanS: 3*NumSubjects
    % HealthyOldmeanS: 3*NumSubjects
    % ...
    % config
    
    
    %dim=5*3
    mean_all = [mean(YoungmeanS, 2)';%avearge over the row;
                  mean(HealthyOldmeanS, 2)';
                  mean(MCIPosmeanS, 2)';
                  mean(MCINegmeanS, 2)';
                  mean(MCIUnkmeanS, 2)'];
    %dim=5*3
    std_all = [std(YoungmeanS,0,2)'./sqrt(size(YoungmeanS,2));
               std(HealthyOldmeanS,0,2)'./sqrt(size(HealthyOldmeanS,2));
               std(MCIPosmeanS,0,2)'./sqrt(size(MCIPosmeanS,2));
               std(MCINegmeanS,0,2)'./sqrt(size(MCINegmeanS,2));
               std(MCIUnkmeanS,0,2)'./sqrt(size(MCIUnkmeanS,2))];

    %% set figure info
    f = figure('visible','off','Position', [100 100 1500 500]);
    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)     
    %%% Color definition %%%
    colorForNochange= config.color_scheme_npg(2,:);
    colorForNoDistalCue= config.color_scheme_npg(3,:);
    colorForNoOpticaFlow= config.color_scheme_npg(4,:);
    colorForAll = [colorForNochange;colorForNoDistalCue;colorForNoOpticaFlow];

    b = bar(mean_all, 'grouped', 'FaceAlpha',0.5, 'LineWidth',1);
    b(1).FaceColor = colorForNochange; %set color for Nochange
    b(2).FaceColor = colorForNoDistalCue; %set color for NoDistalCue
    b(3).FaceColor = colorForNoOpticaFlow; %set color for NoOpticaFlow

    hold on
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(mean_all);
    % Get the x coordinate of the bars
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end  
    
    scatter_facealpha=0.2;
    scatter_edgealpha=0.8;
    %plot scatter for Young
    for i=1:nbars
        scatter(x(i,1), YoungmeanS(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha);
    end

    hold on
    %plot scatter for healthyOld
    for i=1:nbars
        scatter(x(i,2), HealthyOldmeanS(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha);
    end   

    hold on
    %plot scatter for MCIPos
    for i=1:nbars
        scatter(x(i,3), MCIPosmeanS(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha);
    end   

    hold on
    %plot scatter for MCIPos
    for i=1:nbars
        scatter(x(i,4), MCINegmeanS(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha);
    end 

    hold on
    %plot scatter for MCIPos
    for i=1:nbars
        scatter(x(i,5), MCIUnkmeanS(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha);
    end     
    
    hold on
    % Plot the errorbars
    errorbar(x',mean_all,std_all,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 14);  

    hold on
    %add y=1 for reference
    yline(1, 'LineStyle','-.', 'LineWidth', 1, 'Color', config.color_scheme_npg(1,:));
    
    %Further post-processing the figure
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , (1:5),... 
        'XTickLabel'  , ["Young", "HealthyOld", "MCIPos", "MCINeg", "MCIUnk"],...
        'LineWidth'   , .5        );
        %'Ytick'       , [0,0.5,1.0,1.5],...
        %'XLim'        , [0.5, 3.5],...
        %'YLim'        , [0, 1.5],...   
        
    legend(b, {'No change' 'No distal cue' 'No optical flow'}, 'Location','northeast', 'NumColumns',3);
    ylabel("Mean walking speed (m/s)");

    %extract pvalue for group, conditino and interaction to show on the figure 
    group_pvalue = dist_anova_tab{2,7};
    condition_pvalue = dist_anova_tab{3,7};
    interaction_pvalue = dist_anova_tab{4,7};
    str = {['Group Pvalue = ',sprintf('%.2g',group_pvalue)],...
           ['Condition Pvalue = ',sprintf('%.2g',condition_pvalue)],...
           ['Interaction Pvalue = ',sprintf('%.2g',interaction_pvalue)]};
    annotation('textbox',[0.2 0.6 0.3 0.3],'String',str,'FitBoxToText','on');
        
    %% save figure
    exportgraphics(f,config.ResultFolder+"/BarMeanSpeed.png",'Resolution',300);
    exportgraphics(f,config.ResultFolder+"/BarMeanSpeed.pdf",'Resolution',300, 'ContentType','vector');

end