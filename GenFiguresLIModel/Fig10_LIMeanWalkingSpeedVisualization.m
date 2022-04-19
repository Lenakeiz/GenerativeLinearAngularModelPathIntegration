%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default'); %for code reproducibility

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
savefolder = pwd + "/Output/";

%% setting the configuration
config.Speed.alpha = 0.9;                                       %Paramanter for running speed calculation
config.Speed.timeOffsetAfterFlagReach = 2;                      %Time to track after flag reached in seconds 
config.Speed.smoothWindow = 10;                                 % tracking rate should be 10Hz so 4 secs window is 40 datapoints
config.Speed.velocityCutoff = 0.2;                              % velocity cutoff to select only the walking part of the reconstructed velocity
config.Speed.timeOffsetForDetectedTemporalWindow = 0.5;         % time in seconds that will push earlier/ the detected rising edge

config.TrialFilter = 0; %merge all conditions
config.UseGlobalSearch = true;

resultfolder = savefolder+"PaperFigs/Fig4";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end


%% Model fitting using Gamma Base Model
%Model parameter gamma, g3, b, sigma, nu. #params=5
config.ModelName = "LIFull";
config.NumParams = 5;

%% Model fitting for YoungControl data
%% calculating tracking path and transoform data
config.Speed.tresholdForBadParticipantL1Recontruction = 1.55;   % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
YoungControls   = CalculateTrackingPath(YoungControls, config);
%transform data
YoungControls = TransformPaths(YoungControls);
YoungmeanS = getMeanWalkingSpeed(YoungControls);

%% Model fitting for HealthyOld data
config.Speed.tresholdForBadParticipantL1Recontruction = 2.0; 
HealthyControls   = CalculateTrackingPath(HealthyControls, config);
%transform data
HealthyOld = TransformPaths(HealthyControls);
HealthyOldmeanS = getMeanWalkingSpeed(HealthyOld);

%% Model fitting for MCIPos
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCIPos   = CalculateTrackingPath(MCIPos, config);
%transform data
MCIPos = TransformPaths(MCIPos);
MCIPosmeanS = getMeanWalkingSpeed(MCIPos);

%% Model fitting for MCINeg
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCINeg   = CalculateTrackingPath(MCINeg, config);
%transform data
MCINeg = TransformPaths(MCINeg);
MCINegmeanS = getMeanWalkingSpeed(MCINeg);

%% Model fitting for MCIUnk
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
Unknown   = CalculateTrackingPath(Unknown, config);
%transform data
MCIUnk = TransformPaths(Unknown);
MCIUnkmeanS = getMeanWalkingSpeed(MCIUnk);

%% Setting colors for using in plots
ColorPattern; 

%% aBarScatter plot of return distance of all groups and conditions
plotBarScatterOfMeanSpeed(YoungmeanS, HealthyOldmeanS, MCIPosmeanS, MCINegmeanS, MCIUnkmeanS)

%% getting Mean Walking Speed of all participants in a group
function meanS = getMeanWalkingSpeed(Dat)
    numSubjs = length(Dat.Reconstructed);
    meanS = [];
    for subj=1:numSubjs
        %length on leg1
        length1 = table2array(Dat{subj}(:,'L1Real'));
        %duration on leg1
        duration1 = table2array(Dat{subj}(:,'T_L1'));
        meanspeed1 = length1./duration1;

        %length on leg2
        length2 = table2array(Dat{subj}(:,'L2Real'));
        %duration on leg2
        duration2 = table2array(Dat{subj}(:,'T_L2'));
        meanspeed2 = length2./duration2;
        
        %mean speed across trials of one participant
        meanspeed = mean([meanspeed1;meanspeed2]);
        
        meanS = [meanS,meanspeed];
    end
end

%% plot Bar scatter of Return distance
function plotBarScatterOfMeanSpeed(YoungmeanS, HealthyOldmeanS, MCIPosmeanS, MCINegmeanS, MCIUnkmeanS, dist_anova_tab, config)
    % NormedDistYoung: 3*NumSubjects
    % NormedDistHealthyOld: 3*NumSubjects
    % ...
    % config
    
    
    %dim=5*3
    mean_all = [mean(NormedDistYoung, 2)';%avearge over the row;
                  mean(NormedDistHealthyOld, 2)';
                  mean(NormedDistMCIPos, 2)';
                  mean(NormedDistMCINeg, 2)';
                  mean(NormedDistMCIUnk, 2)'];
    %dim=5*3
    std_all = [std(NormedDistYoung,0,2)'./sqrt(size(NormedDistYoung,2));
               std(NormedDistHealthyOld,0,2)'./sqrt(size(NormedDistHealthyOld,2));
               std(NormedDistMCIPos,0,2)'./sqrt(size(NormedDistMCIPos,2));
               std(NormedDistMCINeg,0,2)'./sqrt(size(NormedDistMCINeg,2));
               std(NormedDistMCIUnk,0,2)'./sqrt(size(NormedDistMCIUnk,2))];

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
        scatter(x(i,1), NormedDistYoung(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha);
    end

    hold on
    %plot scatter for healthyOld
    for i=1:nbars
        scatter(x(i,2), NormedDistHealthyOld(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha);
    end   

    hold on
    %plot scatter for MCIPos
    for i=1:nbars
        scatter(x(i,3), NormedDistMCIPos(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha);
    end   

    hold on
    %plot scatter for MCIPos
    for i=1:nbars
        scatter(x(i,4), NormedDistMCINeg(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha);
    end 

    hold on
    %plot scatter for MCIPos
    for i=1:nbars
        scatter(x(i,5), NormedDistMCIUnk(i,:),30, ...
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
    ylabel("Actual distance / correct distance");

    %extract pvalue for group, conditino and interaction to show on the figure 
    group_pvalue = dist_anova_tab{2,7};
    condition_pvalue = dist_anova_tab{3,7};
    interaction_pvalue = dist_anova_tab{4,7};
    str = {['Group Pvalue = ',sprintf('%.2g',group_pvalue)],...
           ['Condition Pvalue = ',sprintf('%.2g',condition_pvalue)],...
           ['Interaction Pvalue = ',sprintf('%.2g',interaction_pvalue)]};
    annotation('textbox',[0.2 0.6 0.3 0.3],'String',str,'FitBoxToText','on');
        
    %% save figure
    exportgraphics(f,config.ResultFolder+"/BarNormedReturnDistance.png",'Resolution',300);
    exportgraphics(f,config.ResultFolder+"/BarNormedReturnDistance.pdf",'Resolution',300, 'ContentType','vector');

end

%% One-way anova on merged MCI
function [anova_tab,multicomp_tab] = OnewayAnovaOnMeanSpeed(YoungmeanS, HealthyOldmeanS, MCIPosmeanS, MCINegmeanS, MCIUnkmeanS, config)

    %load configurations necessary for the script
    resultfolder = config.ResultFolder;
    
    %create storing folder for trajectory if not exist
    savefoldername = resultfolder+"/OnewayAnova/";
    if ~exist(savefoldername, 'dir')
       mkdir(savefoldername);
    end
    
    param_names = ["beta", "bG3", "g2", "g3", 'b', "sigma", "nu"];
    param_nums = length(param_names);
    
    anova_tab = cell(0);
    multicomp_tab = cell(0);

    %for each parameter, calculate the anova and multiple test
    for param_idx=1:param_nums
        
        param_name = param_names(param_idx);
        
        %processing the data into a long numeric vector 
        Params = [MCIPosParams(:,param_idx)',...
                  MCINegParams(:,param_idx)',...
                  MCIUnkParams(:,param_idx)',...
                  HealthyOldParams(:,param_idx)',...
                  YoungParams(:,param_idx)'];
        GroupNames = [string(repmat({'MCIPos'},1,size(MCIPosParams,1))),...
                      string(repmat({'MCINeg'},1,size(MCINegParams,1))),...
                      string(repmat({'MCIUnk'},1,size(MCIUnkParams,1))),...
                      string(repmat({'HealthyOld'},1,size(HealthyOldParams,1))),...
                      string(repmat({'Young'},1,size(YoungParams,1)))];
    
        %Do one-way anova with unbalanced design
        [p,tb1, stats]= anova1(Params, GroupNames, 'display','on');
        anova_tab{param_idx} = tb1;
    
        %Do multiple comparisons on main effect 1
        result = multcompare(stats,'CType','bonferroni');
        multicomp_tab{param_idx} = result;
        title("Multiple comparisons with bonferroni correction of : "+param_name);
        saveas(gcf,savefoldername+"MultiComp"+param_name+".png");
        close(gcf);
    
    end
    
end
