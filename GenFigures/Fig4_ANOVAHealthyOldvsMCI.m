%% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
savefolder = pwd + "/Output/";

%% setting the configuration
config.UseGlobalSearch = true;
resultfolder = savefolder+"PaperFigs/Fig4";
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

%% merge MCI together
AllMCIParams = MergeMCI(AllMCIPosParams, AllMCINegParams, AllMCIUnkParams);

%% TwowayAnova
[anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnovaOn_Young_HealthyOld_MergeMCI(AllYoungParams, AllHealthyOldParams, AllMCIParams, config);

%% BarScatter Plot between HealthyOld and MergedMCI for all Fitted Params
%plotBarScatterOfFittedParam(AllHealthyOldParams, AllMCIParams, config);
plotBoxOfFittedParam(AllHealthyOldParams, AllMCIParams, anova_tab, config);
plotBoxOfFittedParamMergeCondition(AllHealthyOldParams, AllMCIParams, multicomp_tab1, config)

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
function plotBarScatterOfFittedParam(AllHealthyOldParams, AllMCIParams, config)
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
            HealthyOldParam = AllHealthyOldParams{TRIAL_FILTER}(:,ParamIndx);
            [HealthyOld_m, HealthyOld_s] = getMeanSem(HealthyOldParam);
            mean_all(TRIAL_FILTER,1) = HealthyOld_m; 
            sem_all(TRIAL_FILTER,1)=HealthyOld_s;
    
            MCIParam = AllMCIParams{TRIAL_FILTER}(:,ParamIndx);
            [MCI_m, MCI_s] = getMeanSem(MCIParam);
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
        colorForHOld = config.color_scheme_npg(5,:);
        colorForMCI = config.color_scheme_npg(2,:);
    
        b = bar(mean_all, 'grouped', 'FaceAlpha',0.5, 'LineWidth',1);
        b(1).FaceColor = colorForHOld; %set color for HealthyOld
        b(2).FaceColor = colorForMCI; %set color for MCI merge
    
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
            HealthyOldParam = AllHealthyOldParams{TRIAL_FILTER}(:,ParamIndx);
            MCIParam = AllMCIParams{TRIAL_FILTER}(:,ParamIndx);    
    
            %plot scatter for HealthyOld
            for i=1:length(HealthyOldParam)
                jitter_value = 0.0*(rand(1)-0.5);
                scatter(x(1,TRIAL_FILTER)+jitter_value, HealthyOldParam(i),30, ...
                "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForHOld, ...
                'MarkerEdgeAlpha', 1, 'MarkerFaceAlpha', 0.3);
            end
    
            hold on
            %plot scatter for MCI
            for i=1:length(MCIParam)
                jitter_value = 0.0*(rand(1)-0.5);
                scatter(x(2,TRIAL_FILTER)+jitter_value, MCIParam(i),30, ...
                "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForMCI, ...
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
        legend(b, {'HealthyOld' 'MCIMerged'}, 'Location','northwest', 'NumColumns',2);
        %xlabel('G_1','Interpreter','tex'); ylabel('G_2','Interpreter','tex');
        
        %% save figure
        exportgraphics(f,config.ResultFolder+"/Bar_"+StoreName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/Bar_"+StoreName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end

%%
function plotBoxOfFittedParam(AllHealthyOldParams, AllMCIParams, anova_tab, config)
    
    numConds = 3; %3 is the condition number
    ParamName = ["Distance gain gamma", "bG_3", "g_2", "Angular gain g_3", "Anglar bias b", "sigma", "nu"];
    StoreName = ["Gamma", "bG_3", "g_2", "g_3", "b", "sigma", "nu"];
    for ParamIndx=1:length(ParamName)

        HealthyOldParamAllConds = []; %dimension are different, so separate from MCIParamAllConds
        MCIParamAllConds = [];
        for TRIAL_FILTER=1:numConds
            %% extract data
            HealthyOldParam = AllHealthyOldParams{TRIAL_FILTER}(:,ParamIndx);
            HealthyOldParamAllConds = [HealthyOldParamAllConds,HealthyOldParam];

            MCIParam = AllMCIParams{TRIAL_FILTER}(:,ParamIndx);
            MCIParamAllConds = [MCIParamAllConds,MCIParam];
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
        colorForHOld = config.color_scheme_npg(5,:);
        colorForMCI = config.color_scheme_npg(2,:);

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

        %% boxplot for each column in HealthOld
        bp1 = boxplot(HealthyOldParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x-shift_value);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %boxplot for each column in MCIMerge
        bp2 = boxplot(MCIParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x+shift_value);
        set(bp2,'linewidth',box_lineWidth);

        %% Coloring each box
        %findobj first getting the three boxes for MCI (from bp2) from right to left
        %then getting the three boxes for HealthyOld (frm bp1) from right to left
        h = findobj(gca,'Tag','Box'); 
        for i = 1:length(h)
            if i<4  %get the MCI box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForMCI,'FaceAlpha',box_color_transparency);
            else %get the HelthyOld box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForHOld,'FaceAlpha',box_color_transparency);
            end
        end

        %% Adjusting median
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% add scatter plot and the mean of HealthyOld
        num_points = size(HealthyOldParamAllConds,1);
        for i=1:size(HealthyOldParamAllConds,2)
            hold on
            x = i*ones(num_points,1)-shift_value+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, HealthyOldParamAllConds(:,i), scatter_markerSize, ...
                    'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForHOld, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
            %add errorbar
            mean_Hold = mean(HealthyOldParamAllConds(:,i));
            sem_Hold = std(HealthyOldParamAllConds(:,i))./sqrt(length(HealthyOldParamAllConds(:,i)));
            errorbar(i-shift_value,mean_Hold,sem_Hold,'k','LineStyle','None', 'LineWidth', 2);    
            hold on
            %add mean point
            scatter(i-shift_value, mean_Hold, 4*scatter_markerSize, 'd',...
                    'filled','MarkerEdgeColor','k', ...
                    'MarkerFaceColor','w', ...
                    'LineWidth',scatter_marker_edgeWidth);
        end

        %% add scatter plot and the mean of MCI
        num_points = size(MCIParamAllConds,1);
        for i=1:size(MCIParamAllConds,2)
            hold on
            x = i*ones(num_points,1)+shift_value+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, MCIParamAllConds(:,i), scatter_markerSize, ...
                    'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForMCI, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
            %add errorbar
            mean_MCI = mean(MCIParamAllConds(:,i));
            sem_MCI = std(MCIParamAllConds(:,i))./sqrt(length(MCIParamAllConds(:,i)));
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
        legend(gca, {'MCIMerged' 'HealthyOld'}, 'Location','northeast', 'NumColumns',2);
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
function plotBoxOfFittedParamMergeCondition(AllHealthyOldParams, AllMCIParams, multicomp_tab1, config)
    
    numConds = 3; %3 is the condition number
    ParamName = ["Distance gain gamma", "bG_3", "g_2", "Angular gain g_3", "Anglar bias b", "sigma", "nu"];
    StoreName = ["Gamma", "bG_3", "g_2", "g_3", "b", "sigma", "nu"];
    for ParamIndx=1:length(ParamName)

        HealthyOldParamAllConds = []; %dimension are different, so separate from MCIParamAllConds
        MCIParamAllConds = [];

        HealthyOldParamAllCondsMergeinColumn = [];
        MCIParamAllCondsMergeinColumn = [];
        for TRIAL_FILTER=1:numConds
            %% extract data
            HealthyOldParam = AllHealthyOldParams{TRIAL_FILTER}(:,ParamIndx);
            HealthyOldParamAllConds = [HealthyOldParamAllConds,HealthyOldParam];
            HealthyOldParamAllCondsMergeinColumn = [HealthyOldParamAllCondsMergeinColumn;HealthyOldParam];

            MCIParam = AllMCIParams{TRIAL_FILTER}(:,ParamIndx);
            MCIParamAllConds = [MCIParamAllConds,MCIParam];
            MCIParamAllCondsMergeinColumn = [MCIParamAllCondsMergeinColumn;MCIParam];
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
        colorForHOld = config.color_scheme_npg(5,:);
        colorForMCI = config.color_scheme_npg(2,:);

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

        %% boxplot for each column in HealthOld
        bp1 = boxplot(HealthyOldParamAllCondsMergeinColumn, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 1);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %boxplot for each column in MCIMerge
        bp2 = boxplot(MCIParamAllCondsMergeinColumn, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 2);
        set(bp2,'linewidth',box_lineWidth);

        %% Coloring each box
        %findobj first getting the box for MCI (from bp2) 
        %then getting the box for HealthyOld (frm bp1)
        h = findobj(gca,'Tag','Box'); 
        %get the MCI box
        patch(get(h(1),'XData'),get(h(1),'YData'),colorForMCI,'FaceAlpha',box_color_transparency);
        %get the HelthyOld box
        patch(get(h(2),'XData'),get(h(2),'YData'),colorForHOld,'FaceAlpha',box_color_transparency);

        %% Adjusting median
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% add scatter plot and the mean of HealthyOld
        condition_mkr = ['o', 'o', 'o'];
        num_points = size(HealthyOldParamAllConds,1);
        for i=1:size(HealthyOldParamAllConds,2)
            hold on
            x = ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, HealthyOldParamAllConds(:,i), scatter_markerSize, ...
                    'filled', ...
                    condition_mkr(i), ... %marker shape
                    'MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForHOld, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
        end

        %add errorbar
        mean_Hold = mean(HealthyOldParamAllCondsMergeinColumn);
        sem_Hold = std(HealthyOldParamAllCondsMergeinColumn)./sqrt(length(HealthyOldParamAllCondsMergeinColumn));
        errorbar(1,mean_Hold,sem_Hold,'k','LineStyle','None', 'LineWidth', 2);    
        hold on
        %add mean point
        scatter(1, mean_Hold, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);

        %% add scatter plot and the mean of MCI
        num_points = size(MCIParamAllConds,1);
        for i=1:size(MCIParamAllConds,2)
            hold on
            x = 2*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, MCIParamAllConds(:,i), scatter_markerSize, ...
                    'filled', ...
                    condition_mkr(i), ... %marker shape
                    'MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForMCI, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
        end
        %add errorbar
        mean_MCI = mean(MCIParamAllCondsMergeinColumn);
        sem_MCI = std(MCIParamAllCondsMergeinColumn)./sqrt(length(MCIParamAllCondsMergeinColumn));
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
            'XTickLabel'  , {'HealthyOld','MCI'},...
            'LineWidth'   , .5        );
            %'Ytick'       , [0,0.5,1.0,1.5],...
            %'YLim'        , [0, 1.5],...   
        ylabel(ParamName(ParamIndx));
        legend(gca, {'MCIMerged', 'HealthyOld'}, 'Location','northeast', 'NumColumns',2);
        %xlabel('G_1','Interpreter','tex'); ylabel('G_2','Interpreter','tex');

        %extract pvalue for multicomparison of Group effect for showing on the figure
        multicomp_result = multicomp_tab1{ParamIndx};
        Pvalue = multicomp_result(1,6); % MCI vs. HealthyOld see Two-way anova for details
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
function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnovaOn_Young_HealthyOld_MergeMCI(AllYoung, AllHealthyOld, AllMCI, config)

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
        [MCIY, MCIGroupNames, MCIConditionNames]=GroupAndRemoveNaN(AllMCI,param_idx,'MCI');
        [HealthyOldY, HealthyOldGroupNames, HealthyOldConditionNames]=GroupAndRemoveNaN(AllHealthyOld,param_idx,'HealthyOld');
        [YoungY, YoungGroupNames, YoungConditionNames]=GroupAndRemoveNaN(AllYoung,param_idx,'Young');
        
        AllY = [MCIY,HealthyOldY,YoungY];
        AllGroupNames = [MCIGroupNames,HealthyOldGroupNames,YoungGroupNames];
        AllConditionNames = [MCIConditionNames,HealthyOldConditionNames,YoungConditionNames];
    
        %Do two-way anova with unbalanced design
        [p,tb1, stats]= anovan(AllY,{AllGroupNames,AllConditionNames},'model','interaction','varnames',{'Groups','Conditions'},'display','on');
        anova_tab{param_idx} = tb1;
    
        %Do multiple comparisons on main effect 1
        result = multcompare(stats,'Dimension',[1],'CType','bonferroni');
        multicomp_tab1{param_idx} = result;
        title("Multiple comparisons with bonferroni correction of : "+param_name);
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
