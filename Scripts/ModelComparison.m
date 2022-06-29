%% Cleaning variables and set intial seed for code reproducibility
clearvars; close all; clc;
rng('default'); 

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/HowettBrain2019_Dataset.mat');

%% setting the configuration
config.Speed.alpha                                      = 0.9;    % Paramanter for running speed calculation
config.Speed.timeOffsetAfterFlagReach                   = 1.5;    % Time to track after flag reached in seconds 
config.Speed.smoothWindow                               = 10;     % tracking rate should be 10Hz so 4 secs window is 40 datapoints
config.Speed.velocityCutoff                             = 0.2;    % velocity cutoff to select only the walking part of the reconstructed velocity
config.Speed.timeOffsetForDetectedTemporalWindow        = 0.4;    % time in seconds that will push earlier/ the detected rising edge
config.UseGlobalSearch                                  = true;
config.TrackedInboundAngularDeltaT                      = 1;
config.includeStand                                     = false;
config.useweber                                         = false;  % only true when use weber law in simple generative models
config.useOoBtrials                                     = true;

resultfolder = pwd+"/Output/ModelFigures/ModelComparison";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% calculating tracking path and transoform data
config.Speed.tresholdForBadParticipantL1Recontruction = 1.55;   % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
YoungControls                   =   TransformPaths(YoungControls);
YoungControls                   =   CalculateTrackingPath(YoungControls, config);
ManuallyScoringYoung;

%% 1, beta, sigma Model
config.ModelName        =   "beta_sigma_nu";
config.ParamName        =   ["beta", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

[~, ~, ~, ~, ~, ~, ~, beta_sigma_IC] = getResultsAllConditions(YoungControls, config);
[beta_sigma_AIC,beta_sigma_BIC, beta_sigma_NLL] = reformatIC(beta_sigma_IC);

%% 2, g2, g3, nu Model
config.ModelName        =   "g2_g3_sigma_nu";
config.ParamName        =   ["g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

[~, ~, ~, ~, ~, ~, ~, g2_g3_nu_IC] = getResultsAllConditions(YoungControls, config);
[g2_g3_nu_AIC,g2_g3_nu_BIC, g2_g3_nu_NLL] = reformatIC(g2_g3_nu_IC);

%% 3, beta, g2, sigma, nu Model
config.ModelName        =   "beta_g2_sigma_nu";
config.ParamName        =   ["beta", "g2", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

[~, ~, ~, ~, ~, ~, ~, beta_g2_sigma_nu_IC] = getResultsAllConditions(YoungControls, config);
[beta_g2_sigma_nu_AIC, beta_g2_sigma_nu_BIC, beta_g2_sigma_nu_NLL] = reformatIC(beta_g2_sigma_nu_IC);

%% 4, beta, g3, sigma, nu Model
config.ModelName        =   "beta_g3_sigma_nu";
config.ParamName        =   ["beta", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

[~, ~, ~, ~, ~, ~, ~, beta_g3_sigma_nu_IC] = getResultsAllConditions(YoungControls, config);
[beta_g3_sigma_nu_AIC,beta_g3_sigma_nu_BIC, beta_g3_sigma_nu_NLL] = reformatIC(beta_g3_sigma_nu_IC);

%% 5, beta, g2, g3, sigma, nu Model
config.ModelName        =   "beta_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

[~, ~, ~, ~, ~, ~, ~, beta_g2_g3_sigma_nu_IC] = getResultsAllConditions(YoungControls, config);
[beta_g2_g3_sigma_nu_AIC, beta_g2_g3_sigma_nu_BIC,beta_g2_g3_sigma_nu_NLL] = reformatIC(beta_g2_g3_sigma_nu_IC);

%% 6, beta, g2, g3, k3, sigma, nu Model
config.ModelName        =   "beta_g2_g3_k3_sigma_nu";
config.ParamName        =   ["beta", "g2", "g3", "k3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

[~, ~, ~, ~, ~, ~, ~, beta_g2_g3_k3_sigma_nu_IC] = getResultsAllConditions(YoungControls, config);
[beta_g2_g3_k3_sigma_nu_AIC, beta_g2_g3_k3_sigma_nu_BIC,beta_g2_g3_k3_sigma_nu_NLL] = reformatIC(beta_g2_g3_k3_sigma_nu_IC);

%% 3, the Gamma Model with both distance and angle error
config.ModelName        = "GammaModel";
config.ParamName        = ["gamma", "bG3", "g2", "g3", 'b', "sigma", "nu"];
config.NumTotalParams   = length(config.ParamName);
config.NumFreeParams    = 4;
[~, ~, ~, ~, GammaIC] = getResultsAllConditions(YoungControls, config);
[GammaAIC, GammaBIC, GammaNLL] = reformatIC(GammaIC);

%% 4, the G1G2 Model with both distance and angle error
config.ModelName        = "G1G2Model";
config.ParamName        = ["bG1", "bG2", "bG3", "g2", "g3", 'b', "sigma", "nu"];
config.NumTotalParams   = length(config.ParamName);
config.NumFreeParams    = 5;
[~, ~, ~, ~, G1G2IC] = getResultsAllConditions(YoungControls, config);
[G1G2AIC, G1G2BIC, G1G2NLL] = reformatIC(G1G2IC);

%%
ModelNames = {'ConstSpeed', 'IntSpeed', 'Gamma', 'G1G2'};

%% Setting colors for using in plots
ColorPattern; 

%% Box Plot of AIC
ICType = "AIC";
All_AIC = [ConstSpeedAIC, IntSpeedAIC, GammaAIC, G1G2AIC];
plotBoxPlot(All_AIC, ModelNames, ICType, config);


%% Box Plot of BIC
ICType = "BIC";
All_BIC = [ConstSpeedBIC, IntSpeedBIC, GammaBIC, G1G2BIC];
plotBoxPlot(All_BIC, ModelNames, ICType, config);


%% Box Plot of AIC
ICType = "NegLogLikelihood";
All_NLL = [ConstSpeedNLL, IntSpeedNLL, GammaNLL, G1G2NLL];
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