%% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
savefolder = pwd + "/Output/";

%% setting the configuration
config.UseGlobalSearch = true;
resultfolder = savefolder+"PaperFigs/Fig7";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% Model fitting using Gamma Base Model
%Model parameter gamma, g3, b, sigma, nu. #params=5
config.ModelName = "BaseModel";
config.NumParams = 5;

%% Model fitting for YoungControl data
YoungControls = TransformPaths(YoungControls);
[~, AllYoungX, AllYoungDX, AllYoungTheta, ~] = getResultsAllConditions(YoungControls, config);
[NormedDistYoung, NormedAngleYoung] = getNormalizedReturnDistanceandAngle(AllYoungX, AllYoungDX, AllYoungTheta);

%% Model fitting for HealthyOld data
HealthyOld = TransformPaths(HealthyControls);
[~, AllHealthyOldX, AllHealthyOldDX, AllHealthyOldTheta, ~] = getResultsAllConditions(HealthyOld, config);
[NormedDistHealthyOld, NormedAngleHealthyOld] = getNormalizedReturnDistanceandAngle(AllHealthyOldX, AllHealthyOldDX, AllHealthyOldTheta);

%% Model fitting for MCIPos
MCIPos = TransformPaths(MCIPos);
[~, AllMCIPosX, AllMCIPosDX, AllMCIPosTheta, ~] = getResultsAllConditions(MCIPos, config);
[NormedDistMCIPos, NormedAngleMCIPos] = getNormalizedReturnDistanceandAngle(AllMCIPosX, AllMCIPosDX, AllMCIPosTheta);

%% Model fitting for MCINeg
MCINeg = TransformPaths(MCINeg);
[~, AllMCINegX, AllMCINegDX, AllMCINegTheta, ~] = getResultsAllConditions(MCINeg, config);
[NormedDistMCINeg, NormedAngleMCINeg] = getNormalizedReturnDistanceandAngle(AllMCINegX, AllMCINegDX, AllMCINegTheta);

%% Model fitting for MCIUnk
MCIUnk = TransformPaths(Unknown);
[~, AllMCIUnkX, AllMCIUnkDX, AllMCIUnkTheta, ~] = getResultsAllConditions(MCIUnk, config);
[NormedDistMCIUnk, NormedAngleMCIUnk] = getNormalizedReturnDistanceandAngle(AllMCIUnkX, AllMCIUnkDX, AllMCIUnkTheta);

%% Setting colors for using in plots
ColorPattern; 

%% Two-way anova on distance of 5 groups and 3 conditions
config.type = "dist";
[dist_anova_tab,dist_multicomp_tab1,dist_multicomp_tab2, dist_multicomp_tab12] = TwowayAnovaOnDistanceOrAngle(NormedDistYoung, NormedDistHealthyOld, NormedDistMCIPos, NormedDistMCINeg, NormedDistMCIUnk, config);

%% Two-way anova on angle of 5 groups and 3 conditions
config.type = "angle";
[angle_anova_tab,angle_multicomp_tab1,angle_multicomp_tab2, angle_multicomp_tab12] = TwowayAnovaOnDistanceOrAngle(NormedAngleYoung, NormedAngleHealthyOld, NormedAngleMCIPos, NormedAngleMCINeg, NormedAngleMCIUnk, config);

%% BarScatter Plot of return distance of all groups and conditions
plotBarScatterOfReturnDistance(NormedDistYoung, NormedDistHealthyOld, NormedDistMCIPos, NormedDistMCINeg, NormedDistMCIUnk, dist_anova_tab, config);
plotBarScatterOfReturnAngle(NormedAngleYoung, NormedAngleHealthyOld, NormedAngleMCIPos, NormedAngleMCINeg, NormedAngleMCIUnk, angle_anova_tab, config)

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

%% get the normalized return distance and angles for all particpants
% each value is a mean value for per participant
function [NormedDist, NormedAngle] = getNormalizedReturnDistanceandAngle(AllX, AllDX, AllTheta)
    numConds = length(AllX);
    NormedDist = [];
    NormedAngle = [];

    for TRIAL_FILTER=1:numConds
        X = AllX{TRIAL_FILTER};
        DX = AllDX{TRIAL_FILTER};
        Theta = AllTheta{TRIAL_FILTER};

        subjectSize = size(X,2);
        NormedDistAllSubjs = zeros(1,subjectSize);
        NormedAngleAllSubjs = zeros(1,subjectSize);

        for subj=1:subjectSize
            subjX = X{subj};
            subjDX = DX{subj};
            subjTheta = Theta{subj};

            sampleSize = size(subjX,2);
            NormedDistPerSubj = zeros(1,sampleSize);
            NormedAnglePerSubj = zeros(1,sampleSize);
            for tr_id=1:sampleSize
                %normalized return distance
                %actual return distance
                l3 = subjDX{tr_id}(3);
                %correct return distance
                p3 = subjX{tr_id}(3,:);
                correct_l3 = norm(p3);

                normalized_dist = l3/correct_l3;
                NormedDistPerSubj(tr_id) = normalized_dist;

                %normalized return angle
                %actual return angle
                theta3 = subjTheta{tr_id}(3);
                %correct return angle 
                p1 = subjX{tr_id}(1,:);
                p2 = subjX{tr_id}(2,:);
                p3 = subjX{tr_id}(3,:);
                vec1 = p3-p2; vec2 = p1-p3;
                correct_theta3=atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
                correct_theta3 = deg2rad(correct_theta3);

                normalized_angle = theta3/correct_theta3;
                NormedAnglePerSubj(tr_id) = normalized_angle;
            end
            NormedDistAllSubjs(subj)=mean(NormedDistPerSubj);
            NormedAngleAllSubjs(subj)=mean(NormedAnglePerSubj);
        end
        NormedDist = [NormedDist;NormedDistAllSubjs]; %dim=3*NumSubjects
        NormedAngle = [NormedAngle;NormedAngleAllSubjs]; %dim=3*NumSubjects
    end
end

%% plot Bar scatter of Return distance
function plotBarScatterOfReturnDistance(NormedDistYoung, NormedDistHealthyOld, NormedDistMCIPos, ...
                                        NormedDistMCINeg, NormedDistMCIUnk, dist_anova_tab, config)
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

%% plot Bar scatter of Return distance
function plotBarScatterOfReturnAngle(NormedAngleYoung, NormedAngleHealthyOld, NormedAngleMCIPos, ...
                                     NormedAngleMCINeg, NormedAngleMCIUnk, angle_anova_tab, config)
    % NormedDistYoung: 3*NumSubjects
    % NormedDistHealthyOld: 3*NumSubjects
    % ...
    % config
    
    
    %dim=5*3
    mean_all = [mean(NormedAngleYoung, 2)';%avearge over the row;
                  mean(NormedAngleHealthyOld, 2)';
                  mean(NormedAngleMCIPos, 2)';
                  mean(NormedAngleMCINeg, 2)';
                  mean(NormedAngleMCIUnk, 2)'];
    %dim=5*3
    std_all = [std(NormedAngleYoung,0,2)'./sqrt(size(NormedAngleYoung,2));
               std(NormedAngleHealthyOld,0,2)'./sqrt(size(NormedAngleHealthyOld,2));
               std(NormedAngleMCIPos,0,2)'./sqrt(size(NormedAngleMCIPos,2));
               std(NormedAngleMCINeg,0,2)'./sqrt(size(NormedAngleMCINeg,2));
               std(NormedAngleMCIUnk,0,2)'./sqrt(size(NormedAngleMCIUnk,2))];

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
        scatter(x(i,1), NormedAngleYoung(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha);
    end

    hold on
    %plot scatter for healthyOld
    for i=1:nbars
        scatter(x(i,2), NormedAngleHealthyOld(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha);
    end   

    hold on
    %plot scatter for MCIPos
    for i=1:nbars
        scatter(x(i,3), NormedAngleMCIPos(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha);
    end   

    hold on
    %plot scatter for MCIPos
    for i=1:nbars
        scatter(x(i,4), NormedAngleMCINeg(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha);
    end 

    hold on
    %plot scatter for MCIPos
    for i=1:nbars
        scatter(x(i,5), NormedAngleMCIUnk(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha);
    end     
    
    hold on
    % Plot the errorbars
    errorbar(x',mean_all,std_all,'k','LineStyle','None', 'LineWidth', 2, 'CapSize',14);  
    
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
    ylabel("Actual angle / correct angle");

    %extract pvalue for group, conditino and interaction to show on the figure 
    group_pvalue = angle_anova_tab{2,7};
    condition_pvalue = angle_anova_tab{3,7};
    interaction_pvalue = angle_anova_tab{4,7};
    str = {['Group Pvalue = ',sprintf('%.2g',group_pvalue)],...
           ['Condition Pvalue = ',sprintf('%.2g',condition_pvalue)],...
           ['Interaction Pvalue = ',sprintf('%.2g',interaction_pvalue)]};
    annotation('textbox',[0.2 0.6 0.3 0.3],'String',str,'FitBoxToText','on');    
        
    %% save figure
    exportgraphics(f,config.ResultFolder+"/BarNormedReturnAngle.png",'Resolution',300);
    exportgraphics(f,config.ResultFolder+"/BarNormedReturnAngle.pdf",'Resolution',300, 'ContentType','vector');

end

%% Two-way anova on all Groups
function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnovaOnDistanceOrAngle(NormedYoung, NormedHealthyOld, NormedMCIPos, NormedMCINeg, NormedMCIUnk, config)

    %load configurations necessary for the script
    resultfolder = config.ResultFolder;
    type = config.type;
    
    %create storing folder for trajectory if not exist
    savefoldername = resultfolder+"/TwowayAnova_"+type+"/";
    if ~exist(savefoldername, 'dir')
       mkdir(savefoldername);
    end

    %processing the data into a long numeric vector 
    [YoungY, YoungGroupNames, YoungConditionNames]=ReGroupData(NormedYoung,'Young');
    [HealthyOldY, HealthyOldGroupNames, HealthyOldConditionNames]=ReGroupData(NormedHealthyOld,'HealthyOld');
    [MCIPosY, MCIPosGroupNames, MCIPosConditionNames]=ReGroupData(NormedMCIPos,'MCIPos');
    [MCINegY, MCINegGroupNames, MCINegConditionNames]=ReGroupData(NormedMCINeg,'MCIPNeg');
    [MCIUnkY, MCIUnkGroupNames, MCIUnkConditionNames]=ReGroupData(NormedMCIUnk,'MCIPUnk');

    AllY = [MCIPosY,MCINegY,MCIUnkY,YoungY,HealthyOldY];
    AllGroupNames = [MCIPosGroupNames,MCINegGroupNames,MCIUnkGroupNames,YoungGroupNames,HealthyOldGroupNames];
    AllConditionNames = [MCIPosConditionNames,MCINegConditionNames,MCIUnkConditionNames,YoungConditionNames,HealthyOldConditionNames];

    %Do two-way anova with unbalanced design
    [p,anova_tab, stats]= anovan(AllY,{AllGroupNames,AllConditionNames},'model','interaction','varnames',{'Groups','Conditions'},'display','on');

    %Do multiple comparisons on main effect 1
    multicomp_tab1 = multcompare(stats,'Dimension',[1],'CType','bonferroni');
    title("Multiple comparisons with bonferroni correction of");
    saveas(gcf,savefoldername+"MultiCompME1.png");
    close(gcf);

    %Do multiple comparisons on main effect 2
    multicomp_tab2 = multcompare(stats,'Dimension',[2],'CType','bonferroni');
    title("Multiple comparisons with bonferroni correction");
    saveas(gcf,savefoldername+"MultiCompME2.png");
    close(gcf);

    %Do multiple comparisons on main effect 1&2
    multicomp_tab12 = multcompare(stats,'Dimension',[1,2],'CType','bonferroni');
    title("Multiple comparisons with bonferroni correction");
    saveas(gcf,savefoldername+"MultiCompME1ME2.png");
    close(gcf);    
end
