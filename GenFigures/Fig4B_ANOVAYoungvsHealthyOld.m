%% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
savefolder = pwd + "/Output/";

%% setting the configuration
config.UseGlobalSearch = true;
resultfolder = savefolder+"PaperFigs/Fig4B";
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

%% BarScatter Plot between Young and HealthyOld for all Fitted Params
plotBoxOfFittedParam(AllYoungParams, AllHealthyOldParams, anova_tab, config);
plotBoxOfFittedParamMergeCondition(AllYoungParams, AllHealthyOldParams, multicomp_tab1, config)

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
function plotBoxOfFittedParam(AllYoungParams, AllHealthyOldParams, anova_tab, config)
    
    numConds = 3; %3 is the condition number
    ParamName = ["Distance gain gamma", "bG_3", "g_2", "Angular gain g_3", "Anglar bias b", "sigma", "nu"];
    StoreName = ["Gamma", "bG_3", "g_2", "g_3", "b", "sigma", "nu"];
    for ParamIndx=1:length(ParamName)

        YoungParamAllConds = [];
        HealthyOldParamAllConds = []; %dimension are different, so separate from YoungParams
        
        for TRIAL_FILTER=1:numConds
            %% extract data
            YoungParam = AllYoungParams{TRIAL_FILTER}(:,ParamIndx);
            YoungParamAllConds = [YoungParamAllConds,YoungParam];

            HealthyOldParam = AllHealthyOldParams{TRIAL_FILTER}(:,ParamIndx);
            HealthyOldParamAllConds = [HealthyOldParamAllConds,HealthyOldParam];
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
        colorForYoung = config.color_scheme_npg(3,:);
        colorForHOld = config.color_scheme_npg(5,:);
        

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

        %% boxplot for each column in Young
        bp1 = boxplot(YoungParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x-shift_value);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %% boxplot for each column in HealthOld
        bp2 = boxplot(HealthyOldParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x+shift_value);
        set(bp2,'linewidth',box_lineWidth);

        %% Coloring each box
        %findobj first getting the three boxes for HealthyOld (from bp2) from right to left
        %then getting the three boxes for Young (frm bp1) from right to left
        h = findobj(gca,'Tag','Box'); 
        for i = 1:length(h)
            if i<4  %get the HelthyOld box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForHOld,'FaceAlpha',box_color_transparency);
            else %get the Young box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForYoung,'FaceAlpha',box_color_transparency);
            end
        end

        %% Adjusting median
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% add scatter plot and the mean of Young
        num_points = size(YoungParamAllConds,1);
        for i=1:size(YoungParamAllConds,2)
            hold on
            x = i*ones(num_points,1)-shift_value+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, YoungParamAllConds(:,i), scatter_markerSize, ...
                    'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForYoung, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
            %add errorbar
            mean_Young = mean(YoungParamAllConds(:,i));
            sem_Young = std(YoungParamAllConds(:,i))./sqrt(length(YoungParamAllConds(:,i)));
            errorbar(i-shift_value,mean_Young,sem_Young,'k','LineStyle','None', 'LineWidth', 2); 
            hold on
            %add mean point
            scatter(i-shift_value, mean_Young, 4*scatter_markerSize, 'd',...
                    'filled','MarkerEdgeColor','k', ...
                    'MarkerFaceColor','w', ...
                    'LineWidth',scatter_marker_edgeWidth);
        end

        %% add scatter plot and the mean of HealthyOld
        num_points = size(HealthyOldParamAllConds,1);
        for i=1:size(HealthyOldParamAllConds,2)
            hold on
            x = i*ones(num_points,1)+shift_value+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, HealthyOldParamAllConds(:,i), scatter_markerSize, ...
                    'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForHOld, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
            %add errorbar
            mean_Hold = mean(HealthyOldParamAllConds(:,i));
            sem_Hold = std(HealthyOldParamAllConds(:,i))./sqrt(length(HealthyOldParamAllConds(:,i)));
            errorbar(i+shift_value,mean_Hold,sem_Hold,'k','LineStyle','None', 'LineWidth', 2);    
            hold on
            %add mean point
            scatter(i+shift_value, mean_Hold, 4*scatter_markerSize, 'd',...
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
        legend(gca, {'Young', 'HealthyOld'}, 'Location','northeast', 'NumColumns',2);
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
function plotBoxOfFittedParamMergeCondition(AllYoungParams, AllHealthyOldParams, multicomp_tab1, config)
    
    numConds = 3; %3 is the condition number
    ParamName = ["Distance gain gamma", "bG_3", "g_2", "Angular gain g_3", "Anglar bias b", "sigma", "nu"];
    StoreName = ["Gamma", "bG_3", "g_2", "g_3", "b", "sigma", "nu"];
    for ParamIndx=1:length(ParamName)

        YoungParamAllConds = [];
        HealthyOldParamAllConds = []; %dimension are different, so separate from MCIParamAllConds
        
        YoungParamAllCondsMergeinColumn = [];
        HealthyOldParamAllCondsMergeinColumn = [];

        for TRIAL_FILTER=1:numConds
            %% extract data
            YoungParam = AllYoungParams{TRIAL_FILTER}(:,ParamIndx);
            YoungParamAllConds = [YoungParamAllConds,YoungParam];
            YoungParamAllCondsMergeinColumn = [YoungParamAllCondsMergeinColumn;YoungParam];
            
            HealthyOldParam = AllHealthyOldParams{TRIAL_FILTER}(:,ParamIndx);
            HealthyOldParamAllConds = [HealthyOldParamAllConds,HealthyOldParam];
            HealthyOldParamAllCondsMergeinColumn = [HealthyOldParamAllCondsMergeinColumn;HealthyOldParam];
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
        colorForYoung = config.color_scheme_npg(3,:);        
        colorForHOld = config.color_scheme_npg(5,:);

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

        %% boxplot for each column in Young
        bp1 = boxplot(YoungParamAllCondsMergeinColumn, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 1);
        set(bp1,'linewidth',box_lineWidth);
        
        hold on
        %% boxplot for each column in HealthOld
        bp2 = boxplot(HealthyOldParamAllCondsMergeinColumn, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 2);
        set(bp2,'linewidth',box_lineWidth);

        %% Coloring each box
        %findobj first getting the box for HealthyOld (from bp2) 
        %then getting the box for Young (frm bp1)
        h = findobj(gca,'Tag','Box'); 
        %get the HelthyOld box
        patch(get(h(1),'XData'),get(h(1),'YData'),colorForHOld,'FaceAlpha',box_color_transparency);
        %get the Young box
        patch(get(h(2),'XData'),get(h(2),'YData'),colorForYoung,'FaceAlpha',box_color_transparency);

        %% Adjusting median
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% add scatter plot and the mean of Young
        condition_mkr = ['o', 'o', 'o'];
        num_points = size(YoungParamAllConds,1);
        for i=1:size(YoungParamAllConds,2)
            hold on
            x = 1*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, YoungParamAllConds(:,i), scatter_markerSize, ...
                    'filled', ...
                    condition_mkr(i), ... %marker shape
                    'MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForYoung, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
        end
        %add errorbar
        mean_Young = mean(YoungParamAllCondsMergeinColumn);
        sem_Young = std(YoungParamAllCondsMergeinColumn)./sqrt(length(YoungParamAllCondsMergeinColumn));
        errorbar(1,mean_Young,sem_Young,'k','LineStyle','None', 'LineWidth', 2); 
        hold on
        %add mean point
        scatter(1, mean_Young, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);    

        %% add scatter plot and the mean of HealthyOld
        num_points = size(HealthyOldParamAllConds,1);
        for i=1:size(HealthyOldParamAllConds,2)
            hold on
            x = 2*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
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
        errorbar(2,mean_Hold,sem_Hold,'k','LineStyle','None', 'LineWidth', 2);    
        hold on
        %add mean point
        scatter(2, mean_Hold, 2*scatter_markerSize, 'd',...
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
        legend(gca, {'HealthyOld', 'Young'}, 'Location','northeast', 'NumColumns',2);
        %xlabel('G_1','Interpreter','tex'); ylabel('G_2','Interpreter','tex');

        %extract pvalue for multicomparison of Group effect for showing on the figure
        multicomp_result = multicomp_tab1{ParamIndx};
        % 1       2             3
        % MCI     HealthyOld    Young       
        Pvalue = multicomp_result(3,6); % Young vs. HealthyOld see Two-way anova for details
        str = {['P = ',sprintf('%.2g',Pvalue)]};
        annotation('textbox',[0.2 0.6 0.3 0.3],'String',str,'FitBoxToText','on');

        %% save figure
        exportgraphics(f,config.ResultFolder+"/ZMergeCondsBox_"+StoreName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/ZMergeCondsBox_"+StoreName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
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

