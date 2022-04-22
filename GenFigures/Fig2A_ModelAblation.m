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
config.UseGlobalSearch                              = true;

resultfolder = savefolder+"PaperFigs/Fig2_ModelAblation";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% calculating tracking path and transoform data
config.Speed.tresholdForBadParticipantL1Recontruction = 1.55;   % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
YoungControls   = CalculateTrackingPath(YoungControls, config);
YoungControls = TransformPaths(YoungControls);  %transform data

%% 1,the allocentric PI model
config.ModelName        = "AlloModel";
config.ParamName        = ["sigma"];
config.UseGlobalSearch  = false; %only one param, convex error func, no need for GlobalSearch
config.useweber         = false; %only true when use weber law in simple generative models
config.NumTotalParams   = length(config.ParamName);
config.NumFreeParams    = 1;
[~, ~, ~, ~, AlloIC] = getResultsAllConditions(YoungControls, config);
[AlloAIC, AlloBIC, AlloNLL] = reformatIC(AlloIC);

%% 2, the egocentrix PI Model
config.ModelName        = "ConstSpeedModel";
config.ParamName        = ["beta", "bG3", "g2", "g3", 'b', "sigma", "nu"];
config.UseGlobalSearch  = true; %only one param, convex error func, no need for GlobalSearch
config.subtype          = "egoNoise"; %choose from 1,egoNoise / 2, onlyDistErr / 3, onlyAngErr_RGb, 
                                                  %4, onlyAngErr_RGmean / 5, DistAngErr_RGb / 6, DistAngErr_RGmean
config.includeStand     = true;
config.useweber         = false; %only true when use weber law in simple generative models
config.NumTotalParams   = length(config.ParamName);
config.NumFreeParams    = 2;
[~, ~, ~, ~, EgoIC] = getResultsAllConditions(YoungControls, config);
[EgoAIC, EgoBIC, EgoNLL] = reformatIC(EgoIC);

%% 3, the ConstSpeed Model with only distance error
config.ModelName        = "ConstSpeedModel";
config.ParamName        = ["beta", "bG3", "g2", "g3", 'b', "sigma", "nu"];
config.UseGlobalSearch  = true; %only one param, convex error func, no need for GlobalSearch
config.subtype          = "onlyDist"; %choose from 1,egoNoise / 2, onlyDist / 3, onlyAng_RGb, 
                                                      %4, onlyAng_RGmean / 5, DistAng_RGb / 6, DistAng_RGmean
config.includeStand     = true;
config.useweber         = false; %only true when use weber law in simple generative models
config.NumTotalParams   = length(config.ParamName);
config.NumFreeParams    = 3;
[~, ~, ~, ~, onlyDistIC] = getResultsAllConditions(YoungControls, config);
[onlyDistAIC, onlyDistBIC, onlyDistNLL] = reformatIC(onlyDistIC);

%% 4, the ConstSpeed Model with only angle error
config.ModelName        = "ConstSpeedModel";
config.ParamName        = ["beta", "bG3", "g2", "g3", 'b', "sigma", "nu"];
config.UseGlobalSearch  = true; %only one param, convex error func, no need for GlobalSearch
config.subtype          = "onlyAng_RGmean"; %choose from 1,egoNoise / 2, onlyDist / 3, onlyAng_RGb, 
                                                      %4, onlyAng_RGmean / 5, DistAng_RGb / 6, DistAng_RGmean
config.includeStand     = true;
config.useweber         = false; %only true when use weber law in simple generative models
config.NumTotalParams   = length(config.ParamName);
config.NumFreeParams    = 3;
[~, ~, ~, ~, onlyAngIC] = getResultsAllConditions(YoungControls, config);
[onlyAngAIC, onlyAngBIC, onlyAngNLL] = reformatIC(onlyAngIC);

%% 5, the ConstSpeed Model with both distance and angle error
config.ModelName        = "ConstSpeedModel";
config.ParamName        = ["beta", "bG3", "g2", "g3", 'b', "sigma", "nu"];
config.UseGlobalSearch  = true; %only one param, convex error func, no need for GlobalSearch
config.subtype          = "DistAng_RGmean"; %choose from 1,egoNoise / 2, onlyDist / 3, onlyAng_RGb, 
                                                      %4, onlyAng_RGmean / 5, DistAng_RGb / 6, DistAng_RGmean
config.includeStand     = true;
config.useweber         = false; %only true when use weber law in simple generative models
config.NumTotalParams   = length(config.ParamName);
config.NumFreeParams    = 4;
[~, ~, ~, ~, DistAngIC] = getResultsAllConditions(YoungControls, config);
[DistAngAIC, DistAngBIC, DistAngNLL] = reformatIC(DistAngIC);

%%
ModelNames = {'Allo', 'Ego', 'Dist', 'Ang', 'DistAng'};

%% Setting colors for using in plots
ColorPattern; 

%% Box Plot of AIC
ICType = "AIC";
All_AIC = [AlloAIC, EgoAIC, onlyDistAIC, onlyAngAIC, DistAngAIC];
plotBoxPlot(All_AIC, ModelNames, ICType, config);


%% Box Plot of BIC
ICType = "BIC";
All_BIC = [AlloBIC, EgoBIC, onlyDistBIC, onlyAngBIC, DistAngBIC];
plotBoxPlot(All_BIC, ModelNames, ICType, config);


%% Box Plot of AIC
ICType = "NegLogLikelihood";
All_NLL = [AlloNLL, EgoNLL, onlyDistNLL, onlyAngNLL, DistAngNLL];
plotBoxPlot(All_NLL, ModelNames, ICType, config);

%% function for Box plot
function plotBoxPlot(data, ModelNames, ICType, config)
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
%     if ICType=='AIC' | ICType=='BIC'
%         ylimup = 75;
%     else
%         ylimup = 35;
%     end    
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , 1:1:100,... 
        'XLim'        , [0.5,length(ModelNames)+0.5],...
        'LineWidth'   , .5        );
    %'YLim'        , [0,ylimup],...

    %% save figure
    exportgraphics(f,config.ResultFolder+"/"+ICType+".png",'Resolution',300);
    exportgraphics(f,config.ResultFolder+"/"+ICType+".pdf",'Resolution',300,'ContentType','vector');
end