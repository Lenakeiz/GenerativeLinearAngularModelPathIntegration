%% Preparing the data
VAM

%% Setting colors for using in plots
ColorPattern;
resultfolder = savefolder + "PaperFigs/Fig1B";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% Getting Information from results:
[YoungControlsPropDist, YoungControlsPropAng]     = getProportionalLinearAndAngularError(YoungControls);
[HealthyControlsPropDist, HealthyControlsPropAng] = getProportionalLinearAndAngularError(HealthyControls);
[MCIUnkPropDist, MCIUnkPropAng]                   = getProportionalLinearAndAngularError(MCIUnk);
[MCINegPropDist, MCINegPropAng]                   = getProportionalLinearAndAngularError(MCINeg);
[MCIPosPropDist, MCIPosPropAng]                   = getProportionalLinearAndAngularError(MCIPos);

%%
config.type = "Distance";
TwowayAnovaOnDistanceOrAngle(YoungControlsPropDist, HealthyControlsPropDist, MCIPosPropDist, MCINegPropDist, MCIUnkPropDist, config);

%% Plotting quantities
plotBarScatterOfReturnDistance(YoungControlsPropDist, HealthyControlsPropDist, MCIPosPropDist, MCINegPropDist, MCIUnkPropDist, [], config);
plotBarScatterOfReturnAngle   (YoungControlsPropAng , HealthyControlsPropAng , MCIPosPropAng , MCINegPropAng , MCIUnkPropAng , [], config);

%% get the normalized return distance and angles for all particpants
% each value is a mean value for per participant per condition
function [PropDist, PropAng] = getProportionalLinearAndAngularError(Group)

    PropDist = [];
    PropAng  = [];

    numConds = length(Group.Results);
    for TRIAL_FILTER = 1:numConds
        
        CondPropDist  = Group.Results{TRIAL_FILTER}.PropDistErr;
        CondPropAng   = Group.Results{TRIAL_FILTER}.PropAngErr;
        subjectSize   = size(CondPropDist,2);

        CondPropDistMeanSubjs  = zeros(1,subjectSize);
        CondPropAngleMeanSubjs = zeros(1,subjectSize);
        
        for subj=1:subjectSize
            CondPropDistMeanSubjs(1,subj)  = mean(cell2mat(CondPropDist{1,subj}),"omitnan");
            CondPropAngleMeanSubjs(1,subj) = mean(cell2mat(CondPropAng{1,subj}),"omitnan");
        end

        PropDist = [PropDist;CondPropDistMeanSubjs];
        PropAng  = [PropAng;CondPropAngleMeanSubjs];
    
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
    mean_all = [mean(NormedDistYoung, 2,"omitnan")';%avearge over the row;
                  mean(NormedDistHealthyOld, 2,"omitnan")';
                  mean(NormedDistMCIPos, 2,"omitnan")';
                  mean(NormedDistMCINeg, 2,"omitnan")';
                  mean(NormedDistMCIUnk, 2,"omitnan")'];
    %dim=5*3
    std_all = [std(NormedDistYoung,0,2,"omitnan")'      ./sqrt(sum(~isnan(NormedDistYoung),2,"omitnan")');
               std(NormedDistHealthyOld,0,2,"omitnan")' ./sqrt(sum(~isnan(NormedDistHealthyOld),2,"omitnan")');
               std(NormedDistMCIPos,0,2,"omitnan")'     ./sqrt(sum(~isnan(NormedDistMCIPos),2,"omitnan")');
               std(NormedDistMCINeg,0,2,"omitnan")'     ./sqrt(sum(~isnan(NormedDistMCINeg),2,"omitnan")');
               std(NormedDistMCIUnk,0,2,"omitnan")'     ./sqrt(sum(~isnan(NormedDistMCIUnk),2,"omitnan")')];

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
%     group_pvalue = dist_anova_tab{2,7};
%     condition_pvalue = dist_anova_tab{3,7};
%     interaction_pvalue = dist_anova_tab{4,7};
%     str = {['Group Pvalue = ',sprintf('%.2g',group_pvalue)],...
%            ['Condition Pvalue = ',sprintf('%.2g',condition_pvalue)],...
%            ['Interaction Pvalue = ',sprintf('%.2g',interaction_pvalue)]};
%     annotation('textbox',[0.2 0.6 0.3 0.3],'String',str,'FitBoxToText','on');
        
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
    mean_all = [mean(NormedAngleYoung, 2,"omitnan")';%avearge over the row;
                  mean(NormedAngleHealthyOld, 2,"omitnan")';
                  mean(NormedAngleMCIPos, 2,"omitnan")';
                  mean(NormedAngleMCINeg, 2,"omitnan")';
                  mean(NormedAngleMCIUnk, 2,"omitnan")'];
    %dim=5*3
    std_all = [std(NormedAngleYoung,0,2,"omitnan")'      ./sqrt(sum(~isnan(NormedAngleYoung),2,"omitnan")');
               std(NormedAngleHealthyOld,0,2,"omitnan")' ./sqrt(sum(~isnan(NormedAngleHealthyOld),2,"omitnan")');
               std(NormedAngleMCIPos,0,2,"omitnan")'     ./sqrt(sum(~isnan(NormedAngleMCIPos),2,"omitnan")');
               std(NormedAngleMCINeg,0,2,"omitnan")'     ./sqrt(sum(~isnan(NormedAngleMCINeg),2,"omitnan")');
               std(NormedAngleMCIUnk,0,2,"omitnan")'     ./sqrt(sum(~isnan(NormedAngleMCIUnk),2,"omitnan")')];


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

%     %extract pvalue for group, conditino and interaction to show on the figure 
%     group_pvalue = angle_anova_tab{2,7};
%     condition_pvalue = angle_anova_tab{3,7};
%     interaction_pvalue = angle_anova_tab{4,7};
%     str = {['Group Pvalue = ',sprintf('%.2g',group_pvalue)],...
%            ['Condition Pvalue = ',sprintf('%.2g',condition_pvalue)],...
%            ['Interaction Pvalue = ',sprintf('%.2g',interaction_pvalue)]};
%     annotation('textbox',[0.2 0.6 0.3 0.3],'String',str,'FitBoxToText','on');    
        
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

%% 
function [Y, GroupNames, ConditionNames]=ReGroupData(Data,groupname)
%Group all the data from five groups and three conditions into a long
%numeric vector for further two-way anova analysis 
%and also output the factor names 

% Removing nans
isnan_idx = ~isnan(Data(1,:));
nochange = Data(1,isnan_idx);
isnan_idx = ~isnan(Data(2,:));
nodistalcues = Data(1,isnan_idx);
isnan_idx = ~isnan(Data(3,:));
noopticflow = Data(1,isnan_idx);
Y = [nochange, nodistalcues, noopticflow];

numSubjsNoChange     = size(nochange,2);
numSubjsNoDistalCues = size(nodistalcues,2);
numSubjsNoOpticFlow  = size(noopticflow,2);

GroupNames = [string(repmat({groupname},1,numSubjsNoChange)),...
              string(repmat({groupname},1,numSubjsNoDistalCues)),...
              string(repmat({groupname},1,numSubjsNoOpticFlow))];

ConditionNames = [string(repmat({'NoChange'},1 ,numSubjsNoChange)),...
                  string(repmat({'NoDistalCue'},1,numSubjsNoDistalCues)),...
                  string(repmat({'NoOpticFlow'},1,numSubjsNoOpticFlow))];  

end