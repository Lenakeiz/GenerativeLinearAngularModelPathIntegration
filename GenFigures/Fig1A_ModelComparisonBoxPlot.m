%% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
YoungControls = TransformPaths(YoungControls);
savefolder = pwd + "/Output/";

%% setting the configuration
resultfolder = savefolder+"PaperFigs/Fig1A";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% 1,the allocentric PI model without weber's law
config.UseGlobalSearch = false; %only one param, convex error func, no need for GlobalSearch
config.ModelName = "Allo";
config.NumParams = 1;
[~, ~, ~, ~, AllYoungIC_Allo] = getResultsAllConditions(YoungControls, config);

%% 2,the allocentric PI model with weber's law
config.UseGlobalSearch = false; %only one param, convex error func, no need for GlobalSearch
config.ModelName = "AlloWeber";
config.NumParams = 1;
[~, ~, ~, ~, AllYoungIC_AlloWeber] = getResultsAllConditions(YoungControls, config);

%% 3,the egocentric PI model
config.UseGlobalSearch = true;
config.ModelName = "Ego";
config.NumParams = 2;
[~, ~, ~, ~, AllYoungIC_Ego] = getResultsAllConditions(YoungControls, config);

%% 4,the egocentric PI model with weber's law
config.UseGlobalSearch = true;
config.ModelName = "EgoWeber";
config.NumParams = 2;
[~, ~, ~, ~, AllYoungIC_EgoWeber] = getResultsAllConditions(YoungControls, config);

%% 5,Our base model
%no consider of the rotation gain factor in leg2, i.e., gamma, G3=1, g2=1, g3, k3, sigma, nu. #params=6
config.UseGlobalSearch = true;
config.ModelName = "BaseModel";
config.NumParams = 5;
[~, ~, ~, ~, AllYoungIC_Base] = getResultsAllConditions(YoungControls, config);

%%
ModelNames = {'Allo', 'AlloWeber', 'Ego', 'EgoWeber', 'Base'};

%% Setting colors for using in plots
ColorPattern; 

%% Box Plot of AIC
ICType = "AIC";

IC_Allo = getRawIC(AllYoungIC_Allo, ICType);
IC_AlloWeber = getRawIC(AllYoungIC_AlloWeber, ICType);
IC_Ego = getRawIC(AllYoungIC_Ego, ICType);
IC_EgoWeber = getRawIC(AllYoungIC_EgoWeber, ICType);
IC_Base = getRawIC(AllYoungIC_Base, ICType);

IC_AllConds = cell(1,3);

IC_Cond1 = [IC_Allo{1}',IC_AlloWeber{1}',IC_Ego{1}', IC_EgoWeber{1}', IC_Base{1}'];
IC_AllConds{1} = IC_Cond1;
plotBoxPlot(IC_Cond1, ModelNames, "NoChange", ICType, config);

IC_Cond2 = [IC_Allo{2}',IC_AlloWeber{2}',IC_Ego{2}', IC_EgoWeber{2}', IC_Base{2}'];
IC_AllConds{2} = IC_Cond2;
plotBoxPlot(IC_Cond2, ModelNames, "NoDistalCue", ICType, config);

IC_Cond3 = [IC_Allo{3}',IC_AlloWeber{3}',IC_Ego{3}', IC_EgoWeber{3}', IC_Base{3}'];
IC_AllConds{3} = IC_Cond3;
plotBoxPlot(IC_Cond3, ModelNames, "NoOpticalFlow", ICType, config);

%plot all conditions in one figure
plotBoxPlotAllConds(IC_AllConds, ModelNames, ICType, config)

% Box Plot of BIC 
ICType = "BIC";

IC_Allo = getRawIC(AllYoungIC_Allo, ICType);
IC_AlloWeber = getRawIC(AllYoungIC_AlloWeber, ICType);
IC_Ego = getRawIC(AllYoungIC_Ego, ICType);
IC_EgoWeber = getRawIC(AllYoungIC_EgoWeber, ICType);
IC_Base = getRawIC(AllYoungIC_Base, ICType);

IC_AllConds = cell(1,3);

IC_Cond1 = [IC_Allo{1}',IC_AlloWeber{1}',IC_Ego{1}', IC_EgoWeber{1}', IC_Base{1}'];
IC_AllConds{1} = IC_Cond1;
plotBoxPlot(IC_Cond1, ModelNames, "NoChange", ICType, config);

IC_Cond2 = [IC_Allo{2}',IC_AlloWeber{2}',IC_Ego{2}', IC_EgoWeber{2}', IC_Base{2}'];
IC_AllConds{2} = IC_Cond2;
plotBoxPlot(IC_Cond2, ModelNames, "NoDistalCue", ICType, config);

IC_Cond3 = [IC_Allo{3}',IC_AlloWeber{3}',IC_Ego{3}', IC_EgoWeber{3}', IC_Base{3}'];
IC_AllConds{3} = IC_Cond3;
plotBoxPlot(IC_Cond3, ModelNames, "NoOpticalFlow", ICType, config);

%plot all conditions in one figure
plotBoxPlotAllConds(IC_AllConds, ModelNames, ICType, config)

% Box Plot of NegLogLikelihood
ICType = "NegLogLikelihood";

IC_Allo = getRawIC(AllYoungIC_Allo, ICType);
IC_AlloWeber = getRawIC(AllYoungIC_AlloWeber, ICType);
IC_Ego = getRawIC(AllYoungIC_Ego, ICType);
IC_EgoWeber = getRawIC(AllYoungIC_EgoWeber, ICType);
IC_Base = getRawIC(AllYoungIC_Base, ICType);

IC_AllConds = cell(1,3);

IC_Cond1 = [IC_Allo{1}',IC_AlloWeber{1}',IC_Ego{1}', IC_EgoWeber{1}', IC_Base{1}'];
IC_AllConds{1} = IC_Cond1;
plotBoxPlot(IC_Cond1, ModelNames, "NoChange", ICType, config);

IC_Cond2 = [IC_Allo{2}',IC_AlloWeber{2}',IC_Ego{2}', IC_EgoWeber{2}', IC_Base{2}'];
IC_AllConds{2} = IC_Cond2;
plotBoxPlot(IC_Cond2, ModelNames, "NoDistalCue", ICType, config);

IC_Cond3 = [IC_Allo{3}',IC_AlloWeber{3}',IC_Ego{3}', IC_EgoWeber{3}', IC_Base{3}'];
IC_AllConds{3} = IC_Cond3;
plotBoxPlot(IC_Cond3, ModelNames, "NoOpticalFlow", ICType, config);

%plot all conditions in one figure
plotBoxPlotAllConds(IC_AllConds, ModelNames, ICType, config)

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

%% function for get the IC data out from the data structure
function IC_AllConds=getRawIC(IC, ICType)
    IC_AllConds = cell(1,3);
    for i=1:3 %three conditions
        if ICType == "AIC"
            icdata = [];
            for j=1:length(IC{i})
                icdata = [icdata,IC{i}{j}.aic];
            end
        elseif ICType == "BIC"
            icdata = [];
            for j=1:length(IC{i})
                icdata = [icdata,IC{i}{j}.bic];
            end
        elseif ICType == "NegLogLikelihood"
            icdata = [];
            for j=1:length(IC{i})
                icdata = [icdata,IC{i}{j}.negll];
            end      
        else
            error("Choose correct IC type!")
        end
        IC_AllConds{i}=icdata;
    end
end

%% function for Box plot
function plotBoxPlot(data, ModelNames, CondType, ICType, config)
    f = figure('visible','off','Position', [100 100 600 400]);
    
    %%%set paramsters
    num_boxplots = size(data,2);
    box_lineWidth = 0.5;
    whisker_value = 1.5;
    box_widths_value = 0.6;
    box_color_transparency = 0.5; %faceAlpha
    median_lineWidth = 2;
    median_color = 'k';
    scatter_jitter_value = 0.1;
    scatter_markerSize=10;
    scatter_marker_edgeColor = 'k';
    scatter_marker_edgeWidth = 0.5;
    scatter_color_transparency = 0.7; %faceAlpha

    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)
   
    box_colors = config.color_scheme_npg(1:num_boxplots,:);
    
    %% main boxplot one box for each column in data
    bp = boxplot(data, 'whisker',whisker_value,'symbol','', ... %symbol ='' making outlier invisible
                'color','k', ...
                'Notch','on',...
                'widths',box_widths_value,...
                'labels', ModelNames);
    
    set(bp,'linewidth',box_lineWidth);

    if ICType=="AIC"
        ylabel('Akaike information criterion (AIC)');
    elseif ICType=="BIC"
        ylabel('Bayesian Inference Criterion (BIC)');
    elseif ICType=="NegLogLikelihood"
        ylabel('Negative Loglikelihood');
    else
        error("Choose correct IC type!");
    end
    
    %% Coloring each box
    h = findobj(gca,'Tag','Box');
    for i = 1:length(h)
        %note that first getting the most right box, so do "length(h)-i+1"
        patch(get(h(i),'XData'),get(h(i),'YData'),box_colors(i,:),'FaceAlpha',box_color_transparency);
    end
    % Sending patch to back of the figure so that median can be drawn on top of it
    set(gca,'children',flipud(get(gca,'children'))) 

    %% Adjusting median
    h=findobj(gca,'tag','Median');
    for i = 1:length(h)
        h(i).LineWidth = median_lineWidth;
        h(i).Color = median_color;
    end

    %% add scatter plot and the mean of the data
    num_points = size(data,1);
    num_boxplots = size(data,2);
    for i=1:size(data,2)
        hold on
        x = i*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, data(:,i), scatter_markerSize, ...
                'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',box_colors(num_boxplots-i+1,:), ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 
        hold on
        scatter(i, mean(data(:,i)), 4*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);
    end

    %% Further post-processing the figure
    if ICType=='AIC' | ICType=='BIC'
        ylimup = 75;
    else
        ylimup = 35;
    end    
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , 1:1:100,... 
        'XLim'        , [0.5,length(ModelNames)+0.5],...
        'YLim'        , [0,ylimup],...
        'LineWidth'   , .5        );

    %% save figure
    exportgraphics(f,config.ResultFolder+"/"+ICType+"_"+CondType+".png",'Resolution',300);
end

%% function for Box plot of all conditions
function plotBoxPlotAllConds(IC_AllConds, ModelNames, ICType, config)

    f = figure('visible','off','Position', [100 100 1600 400]);
    conditionName = {'No Change', 'No Distal Cue', 'No Optical Flow'};

    nrows = 1;ncols = 3; 
    % w and h of each axis in normalized units
    axisw = (1 / ncols) * 0.9; axish = (1 / nrows) * 0.8;

    for kk=1:3
        % calculate the left, bottom coordinate of this subplot
        axisl = (axisw+0.05) * (kk-1);
        subplot(1,3,kk, 'position', [axisl, 0.1, axisw, axish]);
        
        data = IC_AllConds{kk};

        %%%set paramsters
        num_boxplots = size(data,2);
        box_lineWidth = 0.5;
        whisker_value = 1.5;
        box_widths_value = 0.6;
        box_color_transparency = 0.5; %faceAlpha
        median_lineWidth = 2;
        median_color = 'k';
        scatter_jitter_value = 0.1;
        scatter_markerSize=10;
        scatter_marker_edgeColor = 'k';
        scatter_marker_edgeWidth = 0.5;
        scatter_color_transparency = 0.7; %faceAlpha
    
        %%% Font type and size setting %%%
        % Using Arial as default because all journals normally require the font to
        % be either Arial or Helvetica
        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)
    
        %%% Color definition %%%
        box_colors = config.color_scheme_npg(1:num_boxplots,:);
        
        %% main boxplot one box for each column in data
        bp = boxplot(data, 'whisker',whisker_value,'symbol','', ... %symbol ='' making outlier invisible
                    'color','k', ...
                    'Notch','on',...
                    'widths',box_widths_value,...
                    'labels', ModelNames);
        
        set(bp,'linewidth',box_lineWidth);
    
        if ICType=="AIC"
            ylabel('Akaike information criterion (AIC)');
        elseif ICType=="BIC"
            ylabel('Bayesian Inference Criterion (BIC)');
        elseif ICType=="NegLogLikelihood"
            ylabel('Negative Loglikelihood');
        else
            error("Choose correct IC type!");
        end
        
        %% Coloring each box
        h = findobj(gca,'Tag','Box');
        for i = 1:length(h)
            %note that first getting the most right box, so do "length(h)-i+1"
            patch(get(h(i),'XData'),get(h(i),'YData'),box_colors(i,:),'FaceAlpha',box_color_transparency);
        end
        % Sending patch to back of the figure so that median can be drawn on top of it
        %set(gca,'children',flipud(get(gca,'children'))) 
    
        %% Adjusting median
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end
    
        %% add scatter plot and the mean of the data
        num_points = size(data,1);
        num_boxplots = size(data,2);
        for i=1:size(data,2)
            hold on
            x = i*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, data(:,i), scatter_markerSize, ...
                    'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',box_colors(num_boxplots-i+1,:), ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
            scatter(i, mean(data(:,i)), 4*scatter_markerSize, 'd',...
                    'filled','MarkerEdgeColor','k', ...
                    'MarkerFaceColor','w', ...
                    'LineWidth',scatter_marker_edgeWidth);
        end
    
        %% Further post-processing the figure
        if ICType=='AIC' | ICType=='BIC'
            ylimup = 75;
        else
            ylimup = 35;
        end
        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.1 .1 .1], ...
            'YColor'      , [.1 .1 .1], ...
            'XTick'       , 1:1:100,... 
            'XLim'        , [0.5,length(ModelNames)+0.5],...
            'YLim'        , [0,ylimup],...
            'LineWidth'   , .5        );
        title(conditionName{kk});
        
    end

    %% save figure
    exportgraphics(f,config.ResultFolder+"/AAll_"+ICType+".png",'Resolution',300);
end

