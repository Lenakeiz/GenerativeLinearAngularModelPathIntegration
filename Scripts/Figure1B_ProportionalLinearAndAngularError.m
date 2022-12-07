%% Preparing the data
VAM_PrepareBaseConfig

%% Preprocessing the data
VAM_PreprocessData

%% Setting the model we are interested in
% Eventually modify config paramteters we are interested in. For example
% for this graph we are not interested in running the model so we will
% force it to not run
config.ModelName        =   "beta_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "g2", "g3", "sigma", "nu"];
config.NumParams = 100; % Set this here to prevent the modelling to run
% Run the model
VAM

%% Model completed, preparing the data for plotting figures
config.ResultFolder = pwd + "/Output/PaperFigs/Fig1B";
%create storing folder for trajectory if not exist
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

%% Getting Information from results:
[YoungControlsPropDist, YoungControlsPropAng]     = getProportionalLinearAndAngularError(YoungControls);
[HealthyControlsPropDist, HealthyControlsPropAng] = getProportionalLinearAndAngularError(HealthyControls);
[MCIUnkPropDist, MCIUnkPropAng]                   = getProportionalLinearAndAngularError(MCIUnk);
[MCINegPropDist, MCINegPropAng]                   = getProportionalLinearAndAngularError(MCINeg);
[MCIPosPropDist, MCIPosPropAng]                   = getProportionalLinearAndAngularError(MCIPos);

%%
plotInfo.type = "ProportionalDistance";
[anova_tab_dist,multicomp_tab1_dist,multicomp_tab2_dist, multicomp_tab12_dist]     = TwowayAnovaAllGroupsData(YoungControlsPropDist, HealthyControlsPropDist, MCIPosPropDist, MCINegPropDist, MCIUnkPropDist, config, plotInfo);

%%
plotInfo.type = "ProportionalAngle";
[anova_tab_angle,multicomp_tab1_angle,multicomp_tab2_angle, multicomp_tab12_angle] = TwowayAnovaAllGroupsData(YoungControlsPropAng, HealthyControlsPropAng, MCIPosPropAng, MCINegPropAng, MCIUnkPropAng, config, plotInfo);

%% Plotting
ColorPattern;

%%
plotInfo.defaultTextSize = 20;
plotInfo.defaultLineSize = 2;
plotInfo.barFaceAlpha = 0.5;
plotInfo.barLineWidth = 1.5;
plotInfo.scatterFaceAlpha = 0.2;
plotInfo.scatterEdgeAlpha = 0.8;
plotInfo.scatterDataSize = 32;
plotInfo.errorBarWidth = 2.0;
plotInfo.type = "ProportionalDistance";
plotInfo.YLabel = "Actual distance / correct distance";
plotInfo.visible = "off";
% These values have been set after checking the output of the two way anova
% This adds only stars on top of the groups
plotInfo.addStar = [nan nan nan nan 1];
plotInfo.starXValues = [1 2 3 4 5];
plotInfo.starYValues = [nan nan nan nan 1.5];
plotInfo.starTextSize = [35 35 35 35 35];
% This adds sigstar bars using adjustablesigstar
plotInfo.addSigmaStarbars = 1;
plotInfo.sigmaStarBars = [4 5; 3 5; 2 5; 1 5; 1 2];
plotInfo.sigmaStarBarsPValues = [0.001; 0.001; 0.001; 0.001; 0.02];
plotInfo.sigmaStarLineWidth = 2.5;
plotInfo.sigmaStarTextSize  = 20;
plotInfo.sigmaBarSeparation = 0.075;
plotInfo.yLim = [0 1.75];
plotInfo.ticksStep = 0.5;

%plotBarScatter(YoungControlsPropDist, HealthyControlsPropDist, MCIPosPropDist, MCINegPropDist, MCIUnkPropDist, anova_tab_dist, multicomp_tab1_dist, config, plotInfo);

plotInfo.yLim = [0 1.75];
plotInfo.ticksStep = 0.5;
plotMergedBarScatter(mean(YoungControlsPropDist,1,'omitnan'), mean(HealthyControlsPropDist,1,'omitnan'), mean(MCIPosPropDist,1,'omitnan'), mean(MCINegPropDist,1,'omitnan'), mean(MCIUnkPropDist,1,'omitnan'), anova_tab_dist, multicomp_tab1_dist, config, plotInfo);

%
plotInfo.defaultTextSize = 20;
plotInfo.defaultLineSize = 2;
plotInfo.barFaceAlpha = 0.5;
plotInfo.barLineWidth = 1.5;
plotInfo.scatterFaceAlpha = 0.2;
plotInfo.scatterEdgeAlpha = 0.8;
plotInfo.scatterDataSize = 32;
plotInfo.errorBarWidth = 2.0;
plotInfo.type = "ProportionalAngle";
plotInfo.YLabel = "Actual angle / correct angle";
plotInfo.visible = "off";
% These values have been set after checking the output of the two way anova
% This adds only stars on top of the groups
plotInfo.addStar = [nan nan nan nan 1];
plotInfo.starXValues = [1 2 3 4 5];
plotInfo.starYValues = [nan nan nan nan 2.1];
plotInfo.starTextSize = [35 35 35 35 35];
% This adds sigstar bars using adjustablesigstar
plotInfo.addSigmaStarbars = 1;
plotInfo.sigmaStarBars = [4 5; 3 5; 2 5; 1 5];
plotInfo.sigmaStarBarsPValues = [0.001; 0.001; 0.001; 0.001];
plotInfo.sigmaStarLineWidth = 2.5;
plotInfo.sigmaStarTextSize  = 20;
plotInfo.sigmaBarSeparation = 0.06;
plotInfo.yLim = [0 2.5];
plotInfo.ticksStep = 0.5;

%plotBarScatter(YoungControlsPropAng , HealthyControlsPropAng , MCIPosPropAng , MCINegPropAng , MCIUnkPropAng , anova_tab_angle, multicomp_tab1_angle, config, plotInfo);

plotInfo.yLim = [0 2.5];
plotMergedBarScatter(mean(YoungControlsPropAng,1,'omitnan'), mean(HealthyControlsPropAng,1,'omitnan'), mean(MCIPosPropAng,1,'omitnan'), mean(MCINegPropAng,1,'omitnan'), mean(MCIUnkPropAng,1,'omitnan'), anova_tab_angle, multicomp_tab1_angle, config, plotInfo);

%%
clear plotInfo

%% get the normalized return distance and angles for all particpants
% each value is a mean value for per participant per condition
function [PropDist, PropAng] = getProportionalLinearAndAngularError(Group)

    PropDist = [];
    PropAng  = [];

    for TRIAL_FILTER = 1:3
        
        CondPropDist  = Group.Results.PropDistErr{TRIAL_FILTER};
        CondPropAng   = Group.Results.PropAngErr{TRIAL_FILTER};
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
function plotBarScatter(YoungData, HealthyOldData, MCIPosData, MCINegData, MCIUnkData, ANOVA_tab, ANOVA_multcomp_group, config, plotInfo)
    % YoungData: 3*NumSubjects
    % HealthyOldData: 3*NumSubjects
    % ...
    % config
       
    %dim=5*3
    mean_all = [mean(YoungData, 2,"omitnan")';%avearge over the row;
                  mean(HealthyOldData, 2,"omitnan")';
                  mean(MCIUnkData, 2,"omitnan")';
                  mean(MCINegData, 2,"omitnan")';
                  mean(MCIPosData, 2,"omitnan")'];
    %dim=5*3
    std_all = [std(YoungData,0,2,"omitnan")'      ./sqrt(sum(~isnan(YoungData),2,"omitnan")');
               std(HealthyOldData,0,2,"omitnan")' ./sqrt(sum(~isnan(HealthyOldData),2,"omitnan")');
               std(MCIUnkData,0,2,"omitnan")'     ./sqrt(sum(~isnan(MCIUnkData),2,"omitnan")');
               std(MCINegData,0,2,"omitnan")'     ./sqrt(sum(~isnan(MCINegData),2,"omitnan")');
               std(MCIPosData,0,2,"omitnan")'     ./sqrt(sum(~isnan(MCIPosData),2,"omitnan")')];

    % set figure info
    f = figure('visible',plotInfo.visible,'Position', [100 100 1500 500]);
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

    b = bar(mean_all, 'grouped', 'FaceAlpha',plotInfo.barFaceAlpha, 'LineWidth',plotInfo.barLineWidth);
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
    
    scatter_facealpha=plotInfo.scatterFaceAlpha;
    scatter_edgealpha=plotInfo.scatterEdgeAlpha;
    scatter_datasize =plotInfo.scatterDataSize;
    %plot scatter for Young
    for i=1:nbars
        scatter(x(i,1), YoungData(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha, SizeData=scatter_datasize);
    end

    hold on
    %plot scatter for healthyOld
    for i=1:nbars
        scatter(x(i,2), HealthyOldData(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha, SizeData=scatter_datasize);
    end   

    hold on
    %plot scatter for MCI Unk
    for i=1:nbars
        scatter(x(i,3), MCIUnkData(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha, SizeData=scatter_datasize);
    end

    hold on
    %plot scatter for MCI Neg
    for i=1:nbars
        scatter(x(i,4), MCINegData(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha, SizeData=scatter_datasize);
    end 

    hold on
    %plot scatter for MCIPos
    for i=1:nbars
        scatter(x(i,5), MCIPosData(i,:),30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForAll(i,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha, SizeData=scatter_datasize);
    end   
    
    hold on
    % Plot the errorbars
    errorbar(x',mean_all,std_all,'k','LineStyle','None', 'LineWidth', plotInfo.errorBarWidth, 'CapSize', 14);  

    hold on
    %add y=1 for reference
    yline(1, 'LineStyle','-.', 'LineWidth', 1.5, 'Color', config.color_scheme_npg(1,:), 'Alpha', 0.8);
    
    %Further post-processing the figure
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , (1:5),... 
        'XTickLabel'  , ["Young", "Elderly", "MCI Unk", "MCI Neg", "MCI Pos"],...
        'LineWidth'   , .5        );
        %'Ytick'       , [0,0.5,1.0,1.5],...
        %'XLim'        , [0.5, 3.5],...
        %'YLim'        , [0, 1.5],...   
        
    lg = legend(b, {'No change' 'No distal cue' 'No optical flow'}, 'Location','northwest', 'NumColumns', 3, 'AutoUpdate','off');

    %%%% extract pvalue for group, condition and interaction to show on the figure 
    %     group_pvalue = ANOVA_tab{2,7};
    %     condition_pvalue = ANOVA_tab{3,7};
    %     interaction_pvalue = ANOVA_tab{4,7};
    %     str = {['Group Pvalue = ',sprintf('%.2g',group_pvalue)],...
    %            ['Condition Pvalue = ',sprintf('%.2g',condition_pvalue)],...
    %            ['Interaction Pvalue = ',sprintf('%.2g',interaction_pvalue)]};
    %     annotation('textbox',[0.2 0.6 0.3 0.3],'String',str,'FitBoxToText','on');

    for i=1:length(plotInfo.addStar)
        if(~isnan(plotInfo.addStar(i)))
            text(plotInfo.starXValues(i),plotInfo.starYValues(i),"*",...
                'HorizontalAlignment','Center',...
                'BackgroundColor','none',...
                "FontSize",plotInfo.starTextSize(i));
        end
    end

    ylabel(plotInfo.YLabel);
    ylim(plotInfo.yLim);
    yticks(0:plotInfo.ticksStep:plotInfo.yLim(2));

    hold off;
    
    ax = gca;
    ax.LineWidth = plotInfo.defaultLineSize;

    % save figure
    exportgraphics(f,config.ResultFolder+"/BarNormedReturn" + plotInfo.type + ".png",'Resolution',300);
    exportgraphics(f,config.ResultFolder+"/BarNormedReturn" + plotInfo.type + ".pdf",'Resolution',300, 'ContentType','vector');

end

%% plot Bar scatter of Return distance
function plotMergedBarScatter(YoungData, HealthyOldData, MCIPosData, MCINegData, MCIUnkData, ANOVA_tab, ANOVA_multcomp_group, config, plotInfo)
    % Data has to be formatted to be one line for each single group. Group
    % plotting is always young, older, mci unk, mci neg, mci pos
    % ColorPattern should have been called before this
    % YoungData: 1*NumSubjects
    % HealthyOldData: 1*NumSubjects
    % ...
    % config
    
    % Recombining the data 
    AllData{1} = YoungData; AllData{2} = HealthyOldData; AllData{3} = MCIUnkData; AllData{4} = MCINegData; AllData{5} = MCIPosData;
    %dim=5*3
    mean_all   = [mean(YoungData, 2,"omitnan")';%avearge over the row;
                  mean(HealthyOldData, 2,"omitnan")';
                  mean(MCIUnkData, 2,"omitnan")';
                  mean(MCINegData, 2,"omitnan")';
                  mean(MCIPosData, 2,"omitnan")'];
    %dim=5*3
    std_all = [std(YoungData,0,2,"omitnan")'      ./sqrt(sum(~isnan(YoungData),2,"omitnan")');
               std(HealthyOldData,0,2,"omitnan")' ./sqrt(sum(~isnan(HealthyOldData),2,"omitnan")');
               std(MCIUnkData,0,2,"omitnan")'     ./sqrt(sum(~isnan(MCIUnkData),2,"omitnan")');
               std(MCINegData,0,2,"omitnan")'     ./sqrt(sum(~isnan(MCINegData),2,"omitnan")');
               std(MCIPosData,0,2,"omitnan")'     ./sqrt(sum(~isnan(MCIPosData),2,"omitnan")')];

    % set figure info
    f = figure('visible',plotInfo.visible,'Position', [100 100 850 500]);
    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',plotInfo.defaultTextSize)
    set(0,'DefaultTextFontSize',plotInfo.defaultTextSize)
    
    %% Plotting the bar
    b = bar(mean_all, 'grouped');
    b.FaceColor = 'flat';
    b.FaceAlpha = plotInfo.barFaceAlpha;
    b.LineWidth = plotInfo.barLineWidth;
    for j=1:length(mean_all)
        b.CData(j,:) = config.color_scheme_group(j,:);
    end

    hold on

    % Get the x coordinate of the bars
    x = b.XEndPoints;
    
    scatter_facealpha=plotInfo.scatterFaceAlpha;
    scatter_edgealpha=plotInfo.scatterEdgeAlpha;
    scatter_datasize=plotInfo.scatterDataSize;

    for j=1:length(mean_all)
        %plot scatter for Young
        scatter(x(j), AllData{j},30, ...
        "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor',config.color_scheme_group(j,:), ...
        'MarkerEdgeAlpha', scatter_edgealpha, 'MarkerFaceAlpha', scatter_facealpha, SizeData=scatter_datasize);
    end
    
    hold on
    % Plot the errorbars
    errorbar(x',mean_all,std_all,'k','LineStyle','None', 'LineWidth', plotInfo.errorBarWidth, 'CapSize', 14);  

    hold on
    %add y=1 for reference
    yline(1, 'LineStyle','-.', 'LineWidth', 1.5, 'Color', config.color_scheme_npg(1,:), 'Alpha', 0.8);
    
    %Further post-processing the figure
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , (1:5),... 
        'XTickLabel'  , ["Young", "Elderly", "MCI Unk", "MCI Neg", "MCI Pos"],...
        'LineWidth'   , .5        );
        %'Ytick'       , [0,0.5,1.0,1.5],...
        %'XLim'        , [0.5, 3.5],...
        %'YLim'        , [0, 1.5],...   

    sigstaroptions.textSize      = plotInfo.sigmaStarTextSize;
    sigstaroptions.lineWidth     = plotInfo.sigmaStarLineWidth;
    sigstaroptions.barSeparation = plotInfo.sigmaBarSeparation;

    if(plotInfo.addSigmaStarbars == 1)
        for j = 1:height(plotInfo.sigmaStarBars)
            adjustablesigstar(plotInfo.sigmaStarBars(j,:),plotInfo.sigmaStarBarsPValues(j),0,sigstaroptions);
        end
    end
        
    ylabel(plotInfo.YLabel);
    ylim(plotInfo.yLim);
    yticks(0:plotInfo.ticksStep:plotInfo.yLim(2));
    
    xtickangle(35)

    ax = gca;
    ax.LineWidth = plotInfo.defaultLineSize;

    hold off;
    
    % save figure
    exportgraphics(f,config.ResultFolder+"/BarMergedNormedReturn" + plotInfo.type + ".png",'Resolution',300);
    exportgraphics(f,config.ResultFolder+"/BarMergedNormedReturn" + plotInfo.type + ".pdf",'Resolution',300, 'ContentType','vector');

end

%% Two-way anova on all Groups
function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnovaAllGroupsData(NormedYoung, NormedHealthyOld, NormedMCIPos, NormedMCINeg, NormedMCIUnk, config, plotInfo)

    % Requires data on the format 3 * numSub where 3 is the 
    %load configurations necessary for the script
    resultfolder = config.ResultFolder;
    type = plotInfo.type;
    
    %create storing folder for trajectory if not exist
    savefoldername = resultfolder+"/TwowayAnova_"+type+"/";
    if ~exist(savefoldername, 'dir')
       mkdir(savefoldername);
    end

    %processing the data into a long numeric vector 
    [YoungY, YoungGroupNames, YoungConditionNames]=ReGroupData(NormedYoung,'Young');
    [HealthyOldY, HealthyOldGroupNames, HealthyOldConditionNames]=ReGroupData(NormedHealthyOld,'HealthyOld');
    [MCIPosY, MCIPosGroupNames, MCIPosConditionNames]=ReGroupData(NormedMCIPos,'MCIPos');
    [MCINegY, MCINegGroupNames, MCINegConditionNames]=ReGroupData(NormedMCINeg,'MCINeg');
    [MCIUnkY, MCIUnkGroupNames, MCIUnkConditionNames]=ReGroupData(NormedMCIUnk,'MCIUnk');

    AllY              = [YoungY,             HealthyOldY,             MCIUnkY,             MCINegY,             MCIPosY];
    AllGroupNames     = [YoungGroupNames,    HealthyOldGroupNames,    MCIUnkGroupNames,    MCINegGroupNames,    MCIPosGroupNames];
    AllConditionNames = [YoungConditionNames,HealthyOldConditionNames,MCIUnkConditionNames,MCINegConditionNames,MCIPosConditionNames];

    %Do two-way anova with unbalanced design
    [p, anova_tab, stats]= anovan(AllY,{AllGroupNames,AllConditionNames},'model','interaction','varnames',{'Groups','Conditions'},'display','on');

    %Do multiple comparisons on main effect 1
    multicomp_tab1 = multcompare(stats,'Dimension', [1], 'CType','bonferroni');
    title("Multiple comparisons with bonferroni correction (Group)");
    saveas(gcf,savefoldername + "MultiCompME1.png");
    close(gcf);

    %Do multiple comparisons on main effect 2
    multicomp_tab2 = multcompare(stats,'Dimension',[2],'CType','bonferroni');
    title("Multiple comparisons with bonferroni correction (Condition)");
    saveas(gcf,savefoldername+"MultiCompME2.png");
    close(gcf);

    %Do multiple comparisons on main effect 1&2
    multicomp_tab12 = multcompare(stats,'Dimension',[1,2],'CType','bonferroni');
    title("Multiple comparisons with bonferroni correction (Group*Cond)");
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