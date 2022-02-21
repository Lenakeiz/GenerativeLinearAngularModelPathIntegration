%% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
savefolder = pwd + "/Output/";

%% setting the configuration
config.UseGlobalSearch = true;
resultfolder = savefolder+"PaperFigs/Fig6";
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
[AllYoungParams, ~, ~, ~, ~] = getResultsAllConditions(YoungControls, config);

%% Model fitting for HealthyOld data
HealthyOld = TransformPaths(HealthyControls);
[AllHealthyOldParams, ~, ~, ~, ~] = getResultsAllConditions(HealthyOld, config);

%% Model fitting for MCIPos
MCIPos = TransformPaths(MCIPos);
[AllMCIPosParams,  ~, ~, ~, ~] = getResultsAllConditions(MCIPos, config);

%% Model fitting for MCINeg
MCINeg = TransformPaths(MCINeg);
[AllMCINegParams,  ~, ~, ~, ~] = getResultsAllConditions(MCINeg, config);

%% Model fitting for MCIUnk
MCIUnk = TransformPaths(Unknown);
[AllMCIUnkParams,  ~, ~, ~, ~] = getResultsAllConditions(MCIUnk, config);

%% Setting colors for using in plots
ColorPattern; 

%% TwowayAnova
[anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnovaOn_allGroups(AllYoungParams, AllHealthyOldParams, AllMCIPosParams, AllMCINegParams, AllMCIUnkParams, config);

%% BarScatter Plot between MCIPos and MCINeg for all Fitted Params
%plotBarScatterOfFittedParam(AllMCIPosParams, AllMCINegParams, config);
plotBoxOfFittedParam(AllMCIPosParams, AllMCINegParams, anova_tab, config);
plotBoxOfFittedParamMergeCondition(AllMCIPosParams, AllMCINegParams, multicomp_tab1, config)

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

%%
function plotBarScatterOfFittedParam(AllMCIPosParams, AllMCINegParams, config)
    % AllParams: estimated parameter values
    % ParamIndx
    % config
    
    numConds = 3; %3 is the condition number
    mean_all = zeros(numConds, 2); %2 is the group number
    sem_all = zeros(numConds, 2); %2 is the group number

    ParamName = ["Distance gain gamma", "bG_3", "g_2", "Angle gain g_3", "Angle bias b", "sigma", "nu"];
    StoreName = ["Gamma", "bG_3", "g_2", "g_3", "b", "sigma", "nu"];

    for ParamIndx=1:length(ParamName)

        for TRIAL_FILTER=1:numConds
            %% extract data
            MCIPosParam = AllMCIPosParams{TRIAL_FILTER}(:,ParamIndx);
            [MCIPos_m, MCINeg_s] = getMeanSem(MCIPosParam);
            mean_all(TRIAL_FILTER,1) = MCIPos_m; 
            sem_all(TRIAL_FILTER,1)=MCINeg_s;
    
            MCINegParam = AllMCINegParams{TRIAL_FILTER}(:,ParamIndx);
            [MCI_m, MCI_s] = getMeanSem(MCINegParam);
            mean_all(TRIAL_FILTER,2) = MCI_m; 
            sem_all(TRIAL_FILTER,2)=MCI_s;       
        end
    
        %% set figure info
        f = figure('visible','off','Position', [100 100 1000 500]);
        %%% Font type and size setting %%%
        % Using Arial as default because all journals normally require the font to
        % be either Arial or Helvetica
        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)     
        %%% Color definition %%%
        colorForMCIPos = config.color_scheme_npg(6,:);
        colorForMCINeg = config.color_scheme_npg(3,:);
    
        b = bar(mean_all, 'grouped', 'FaceAlpha',0.5, 'LineWidth',1);
        b(1).FaceColor = colorForMCIPos; %set color for MCIPos
        b(2).FaceColor = colorForMCINeg; %set color for MCINeg
    
        hold on
        % Calculate the number of groups and number of bars in each group
        [ngroups,nbars] = size(mean_all);
        % Get the x coordinate of the bars
        x = nan(nbars, ngroups);
        for i = 1:nbars
            x(i,:) = b(i).XEndPoints;
        end
    
        % scatter
        for TRIAL_FILTER=1:numConds
            %% extract data
            MCIPosParam = AllMCIPosParams{TRIAL_FILTER}(:,ParamIndx);
            MCINegParam = AllMCINegParams{TRIAL_FILTER}(:,ParamIndx);    
    
            %plot scatter for MCIPos
            for i=1:length(MCIPosParam)
                jitter_value = 0.0*(rand(1)-0.5);
                scatter(x(1,TRIAL_FILTER)+jitter_value, MCIPosParam(i),30, ...
                "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForMCIPos, ...
                'MarkerEdgeAlpha', 1, 'MarkerFaceAlpha', 0.3);
            end
    
            hold on
            %plot scatter for MCI
            for i=1:length(MCINegParam)
                jitter_value = 0.0*(rand(1)-0.5);
                scatter(x(2,TRIAL_FILTER)+jitter_value, MCINegParam(i),30, ...
                "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForMCINeg, ...
                'MarkerEdgeAlpha', 1, 'MarkerFaceAlpha', 0.3);  
            end
    
        end    
    
        hold on
        % Plot the errorbars
        errorbar(x',mean_all,sem_all,'k','LineStyle','None', 'LineWidth', 2);    
    
        %Further post-processing the figure
        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.1 .1 .1], ...
            'YColor'      , [.1 .1 .1], ...
            'XTick'       , [1:3],... 
            'XTickLabel'  , {'No Change','No Distal Cue', 'No Optical Flow'},...
            'LineWidth'   , .5        );
            %'Ytick'       , [0,0.5,1.0,1.5],...
            %'XLim'        , [0.5, 3.5],...
            %'YLim'        , [0, 1.5],...   
        
        ylabel(ParamName(ParamIndx));
        %legend(b, {'MCIPos' 'MCINeg'}, 'Location','northwest', 'NumColumns',2);
        %xlabel('G_1','Interpreter','tex'); ylabel('G_2','Interpreter','tex');
        
        %% save figure
        exportgraphics(f,config.ResultFolder+"/Bar_"+StoreName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/Bar_"+StoreName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end

%%
function plotBoxOfFittedParam(AllMCIPosParams, AllMCINegParams, anova_tab, config)
    
    numConds = 3; %3 is the condition number
    ParamName = ["Distance gain gamma", "bG_3", "g_2", "Angular gain g_3", "Anglar bias b", "sigma", "nu"];
    StoreName = ["Gamma", "bG_3", "g_2", "g_3", "b", "sigma", "nu"];
    for ParamIndx=1:length(ParamName)

        MCIPosParamAllConds = []; %dimension are different, so separate from MCIParamAllConds
        MCINegParamAllConds = [];
        for TRIAL_FILTER=1:numConds
            %% extract data
            MCIPosParam = AllMCIPosParams{TRIAL_FILTER}(:,ParamIndx);
            MCIPosParamAllConds = [MCIPosParamAllConds,MCIPosParam];

            MCINegParam = AllMCINegParams{TRIAL_FILTER}(:,ParamIndx);
            MCINegParamAllConds = [MCINegParamAllConds,MCINegParam];
        end
    
        %% set figure info
        f = figure('visible','off','Position', [100 100 1000 500]);
        %%% Font type and size setting %%%
        % Using Arial as default because all journals normally require the font to
        % be either Arial or Helvetica
        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)     
        %%% Color definition %%%
        colorForMCIPos = config.color_scheme_npg(6,:);
        colorForMCINeg = config.color_scheme_npg(3,:);

        %set params
        whisker_value = 1.5;
        box_lineWidth = 0.3;
        box_widths_value = 0.2;
        box_color_transparency = 0.5; %faceAlpha
        %center of box (three conditions)
        center_x = [1,2,3];
        shift_value = 0.2; %box shift from center
        median_lineWidth = 2;
        median_color = 'k';
        scatter_jitter_value = 0.1;
        scatter_markerSize=10;
        scatter_marker_edgeColor = 'k';
        scatter_marker_edgeWidth = 0.5;
        scatter_color_transparency = 0.7; %faceAlpha        

        %% boxplot for each column in MCIPos
        bp1 = boxplot(MCIPosParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x-shift_value);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %boxplot for each column in MCIMerge
        bp2 = boxplot(MCINegParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x+shift_value);
        set(bp2,'linewidth',box_lineWidth);

        %% Coloring each box
        %findobj first getting the three boxes for MCINeg (from bp2) from right to left
        %then getting the three boxes for MCIPos (frm bp1) from right to left
        h = findobj(gca,'Tag','Box'); 
        for i = 1:length(h)
            if i<4  %get the MCI box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForMCINeg,'FaceAlpha',box_color_transparency);
            else %get the HelthyOld box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForMCIPos,'FaceAlpha',box_color_transparency);
            end
        end

        %% Adjusting median
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% add scatter plot and the mean of MCIPos
        num_points = size(MCIPosParamAllConds,1);
        for i=1:size(MCIPosParamAllConds,2)
            hold on
            x = i*ones(num_points,1)-shift_value+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, MCIPosParamAllConds(:,i), scatter_markerSize, ...
                    'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForMCIPos, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
            %add errorbar
            mean_MCIPos = mean(MCIPosParamAllConds(:,i));
            sem_MCIPos = std(MCIPosParamAllConds(:,i))./sqrt(length(MCIPosParamAllConds(:,i)));
            errorbar(i-shift_value,mean_MCIPos,sem_MCIPos,'k','LineStyle','None', 'LineWidth', 2);    
            hold on
            %add mean point
            scatter(i-shift_value, mean_MCIPos, 4*scatter_markerSize, 'd',...
                    'filled','MarkerEdgeColor','k', ...
                    'MarkerFaceColor','w', ...
                    'LineWidth',scatter_marker_edgeWidth);
        end

        %% add scatter plot and the mean of MCINeg
        num_points = size(MCINegParamAllConds,1);
        for i=1:size(MCINegParamAllConds,2)
            hold on
            x = i*ones(num_points,1)+shift_value+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, MCINegParamAllConds(:,i), scatter_markerSize, ...
                    'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForMCINeg, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
            %add errorbar
            mean_MCI = mean(MCINegParamAllConds(:,i));
            sem_MCI = std(MCINegParamAllConds(:,i))./sqrt(length(MCINegParamAllConds(:,i)));
            errorbar(i+shift_value,mean_MCI,sem_MCI,'k','LineStyle','None', 'LineWidth', 2); 
            hold on
            %add mean point
            scatter(i+shift_value, mean_MCI, 4*scatter_markerSize, 'd',...
                    'filled','MarkerEdgeColor','k', ...
                    'MarkerFaceColor','w', ...
                    'LineWidth',scatter_marker_edgeWidth);
        end

        %% Further post-processing the figure
        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.1 .1 .1], ...
            'YColor'      , [.1 .1 .1], ...
            'XTick'       , (1:3),... 
            'XLim'        , [0.5, 3.5],...
            'XTickLabel'  , {'No Change','No Distal Cue', 'No Optical Flow'},...
            'LineWidth'   , .5        );
            %'Ytick'       , [0,0.5,1.0,1.5],...
            %'YLim'        , [0, 1.5],...   
        ylabel(ParamName(ParamIndx));
        legend(gca, {'MCIPos' 'MCINeg'}, 'Location','northeast', 'NumColumns',2);
        %xlabel('G_1','Interpreter','tex'); ylabel('G_2','Interpreter','tex');

        %extract pvalue for group, conditino and interaction to show on the figure 
        anova_result = anova_tab{ParamIndx};
        group_pvalue = anova_result{2,7};
        condition_pvalue = anova_result{3,7};
        interaction_pvalue = anova_result{4,7};
        str = {['Group Pvalue = ',sprintf('%.2g',group_pvalue)],...
               ['Condition Pvalue = ',sprintf('%.2g',condition_pvalue)],...
               ['Interaction Pvalue = ',sprintf('%.2g',interaction_pvalue)]};
        annotation('textbox',[0.2 0.6 0.3 0.3],'String',str,'FitBoxToText','on');

        %% save figure
        exportgraphics(f,config.ResultFolder+"/Box_"+StoreName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/Box_"+StoreName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end

%%
function plotBoxOfFittedParamMergeCondition(AllMCIPosParams, AllMCINegParams, multicomp_tab1, config)
    
    numConds = 3; %3 is the condition number
    ParamName = ["Distance gain gamma", "bG_3", "g_2", "Angular gain g_3", "Anglar bias b", "sigma", "nu"];
    StoreName = ["Gamma", "bG_3", "g_2", "g_3", "b", "sigma", "nu"];
    for ParamIndx=1:length(ParamName)

        MCIPosParamAllConds = []; %dimension are different, so separate from MCIParamAllConds
        MCINegParamAllConds = [];

        MCIPosParamAllCondsMergeinColumn = [];
        MCINegParamAllCondsMergeinColumn = [];
        for TRIAL_FILTER=1:numConds
            %% extract data
            MCIPosParam = AllMCIPosParams{TRIAL_FILTER}(:,ParamIndx);
            MCIPosParamAllConds = [MCIPosParamAllConds,MCIPosParam];
            MCIPosParamAllCondsMergeinColumn = [MCIPosParamAllCondsMergeinColumn;MCIPosParam];

            MCINegParam = AllMCINegParams{TRIAL_FILTER}(:,ParamIndx);
            MCINegParamAllConds = [MCINegParamAllConds,MCINegParam];
            MCINegParamAllCondsMergeinColumn = [MCINegParamAllCondsMergeinColumn;MCINegParam];
        end
    
        %% set figure info
        f = figure('visible','off','Position', [100 100 500 500]);
        %%% Font type and size setting %%%
        % Using Arial as default because all journals normally require the font to
        % be either Arial or Helvetica
        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)     
        %%% Color definition %%%
        colorForMCIPos = config.color_scheme_npg(6,:);
        colorForMCINeg = config.color_scheme_npg(3,:);

        %set params
        whisker_value = 1.5;
        box_lineWidth = 0.3;
        box_widths_value = 0.2;
        box_color_transparency = 0.5; %faceAlpha
        median_lineWidth = 2;
        median_color = 'k';
        scatter_jitter_value = 0.1;
        scatter_markerSize=10;
        scatter_marker_edgeColor = 'k';
        scatter_marker_edgeWidth = 0.5;
        scatter_color_transparency = 0.7; %faceAlpha        

        %% boxplot for each column in MCIPOs
        bp1 = boxplot(MCIPosParamAllCondsMergeinColumn, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 1);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %boxplot for each column in MCINeg
        bp2 = boxplot(MCINegParamAllCondsMergeinColumn, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 2);
        set(bp2,'linewidth',box_lineWidth);

        %% Coloring each box
        %findobj first getting the box for MCINeg (from bp2) 
        %then getting the box for MCIPos(frm bp1)
        h = findobj(gca,'Tag','Box'); 
        %get the MCI box
        patch(get(h(1),'XData'),get(h(1),'YData'),colorForMCINeg,'FaceAlpha',box_color_transparency);
        %get the HelthyOld box
        patch(get(h(2),'XData'),get(h(2),'YData'),colorForMCIPos,'FaceAlpha',box_color_transparency);

        %% Adjusting median
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% add scatter plot and the mean of MCIPos
        condition_mkr = ['o', 'o', 'o'];
        num_points = size(MCIPosParamAllConds,1);
        for i=1:size(MCIPosParamAllConds,2)
            hold on
            x = ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, MCIPosParamAllConds(:,i), scatter_markerSize, ...
                    'filled', ...
                    condition_mkr(i), ... %marker shape
                    'MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForMCIPos, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
        end

        %add errorbar
        mean_Hold = mean(MCIPosParamAllCondsMergeinColumn);
        sem_Hold = std(MCIPosParamAllCondsMergeinColumn)./sqrt(length(MCIPosParamAllCondsMergeinColumn));
        errorbar(1,mean_Hold,sem_Hold,'k','LineStyle','None', 'LineWidth', 2);    
        hold on
        %add mean point
        scatter(1, mean_Hold, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);

        %% add scatter plot and the mean of MCI
        num_points = size(MCINegParamAllConds,1);
        for i=1:size(MCINegParamAllConds,2)
            hold on
            x = 2*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, MCINegParamAllConds(:,i), scatter_markerSize, ...
                    'filled', ...
                    condition_mkr(i), ... %marker shape
                    'MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForMCINeg, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
        end
        %add errorbar
        mean_MCI = mean(MCINegParamAllCondsMergeinColumn);
        sem_MCI = std(MCINegParamAllCondsMergeinColumn)./sqrt(length(MCINegParamAllCondsMergeinColumn));
        errorbar(2,mean_MCI,sem_MCI,'k','LineStyle','None', 'LineWidth', 2); 
        hold on
        %add mean point
        scatter(2, mean_MCI, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);        

        %% Further post-processing the figure
        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.1 .1 .1], ...
            'YColor'      , [.1 .1 .1], ...
            'XTick'       , (1:2),... 
            'XLim'        , [0.5, 2.5],...
            'XTickLabel'  , {'MCIPos','MCINeg'},...
            'LineWidth'   , .5        );
            %'Ytick'       , [0,0.5,1.0,1.5],...
            %'YLim'        , [0, 1.5],...   
        ylabel(ParamName(ParamIndx));
        legend(gca, {'MCINeg', 'MCIPos'}, 'Location','northeast', 'NumColumns',2);
        %xlabel('G_1','Interpreter','tex'); ylabel('G_2','Interpreter','tex');

        %extract pvalue for multicomparison of Group effect for showing on the figure
        multicomp_result = multicomp_tab1{ParamIndx};
        Pvalue = multicomp_result(1,6); % MCIPos vs. MCINeg see Two-way anova for details
        str = {['P = ',sprintf('%.2g',Pvalue)]};
        annotation('textbox',[0.2 0.6 0.3 0.3],'String',str,'FitBoxToText','on');

        %% save figure
        exportgraphics(f,config.ResultFolder+"/ZMergeCondsBox_"+StoreName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/ZMergeCondsBox_"+StoreName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end

%%
function [m, s] = getMeanSem(Param)
    m = mean(Param);
    s = std(Param)./sqrt(length(Param));
end

%% Two-way anova on merged MCI
function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnovaOn_allGroups(AllYoung, AllHealthyOld, AllMCIPos, AllMCINeg, AllMCIUnk, config)

    %load configurations necessary for the script
    resultfolder = config.ResultFolder;
    
    %create storing folder for trajectory if not exist
    savefoldername = resultfolder+"/TwowayAnova/";
    if ~exist(savefoldername, 'dir')
       mkdir(savefoldername);
    end
    
    param_names = ["gamma", "bG3", "g2", "g3", 'b', "sigma", "nu"];
    param_nums = length(param_names);
    
    anova_tab = cell(0);
    multicomp_tab1 = cell(0);
    multicomp_tab2 = cell(0);
    multicomp_tab12 = cell(0);
    %for each parameter, calculate the anova and multiple test
    for param_idx=1:param_nums
        
        param_name = param_names(param_idx);
        
    %processing the data into a long numeric vector 
    [YoungY, YoungGroupNames, YoungConditionNames]=GroupAndRemoveNaN(AllYoung,param_idx,'Young');
    [HealthyOldY, HealthyOldGroupNames, HealthyOldConditionNames]=GroupAndRemoveNaN(AllHealthyOld,param_idx,'HealthyOld');
    [MCIPosY, MCIPosGroupNames, MCIPosConditionNames]=GroupAndRemoveNaN(AllMCIPos,param_idx,'MCIPos');
    [MCINegY, MCINegGroupNames, MCINegConditionNames]=GroupAndRemoveNaN(AllMCINeg,param_idx,'MCIPNeg');
    [MCIUnkY, MCIUnkGroupNames, MCIUnkConditionNames]=GroupAndRemoveNaN(AllMCIUnk,param_idx,'MCIPUnk');

    AllY = [MCIPosY,MCINegY,MCIUnkY,YoungY,HealthyOldY];
    AllGroupNames = [MCIPosGroupNames,MCINegGroupNames,MCIUnkGroupNames,YoungGroupNames,HealthyOldGroupNames];
    AllConditionNames = [MCIPosConditionNames,MCINegConditionNames,MCIUnkConditionNames,YoungConditionNames,HealthyOldConditionNames];

    %Do two-way anova with unbalanced design
    [p,tb1, stats]= anovan(AllY,{AllGroupNames,AllConditionNames},'model','interaction','varnames',{'Groups','Conditions'},'display','on');
    anova_tab{param_idx} = tb1;

    %Do multiple comparisons on main effect 1
    result = multcompare(stats,'Dimension',[1],'CType','bonferroni');
    multicomp_tab1{param_idx} = result;
    title("Multiple comparisons with bonferroni correction of parameter: "+param_name);
    saveas(gcf,savefoldername+"MultiCompME1_"+param_name+".png");
    close(gcf);

    %Do multiple comparisons on main effect 2
    result = multcompare(stats,'Dimension',[2],'CType','bonferroni');
    multicomp_tab2{param_idx} = result;
    title("Multiple comparisons with bonferroni correction of parameter: "+param_name);
    saveas(gcf,savefoldername+"MultiCompME2_"+param_name+".png");
    close(gcf);

    %Do multiple comparisons on main effect 1&2
    result = multcompare(stats,'Dimension',[1,2],'CType','bonferroni');
    multicomp_tab12{param_idx} = result;
    title("Multiple comparisons with bonferroni correction of parameter: "+param_name);
    saveas(gcf,savefoldername+"MultiCompME1ME2_"+param_name+".png");
    close(gcf);    
    end
    
end
