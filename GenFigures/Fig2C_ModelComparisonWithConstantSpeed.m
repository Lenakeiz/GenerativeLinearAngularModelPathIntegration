% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
savefolder = pwd + "/Output/";

%% set up the configuration
config.Speed.alpha = 0.9;                                       %Paramanter for running speed calculation
config.Speed.timeOffsetAfterFlagReach = 2;                      %Time to track after flag reached in seconds 
config.Speed.smoothWindow = 10;                                 % tracking rate should be 10Hz so 4 secs window is 40 datapoints
config.Speed.velocityCutoff = 0.1;                              % velocity cutoff to select only the walking part of the reconstructed velocity
config.Speed.timeOffsetForDetectedTemporalWindow = 1.0;         % time in seconds that will push earlier/ the detected rising edge

config.TrialFilter = 0; %merge all conditions
config.UseGlobalSearch = true;

resultfolder = savefolder+"PaperFigs/Fig1C";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

% calculating tracking path and transoform data
config.Speed.tresholdForBadParticipantL1Recontruction = 1.55;   % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
YoungControls   = CalculateTrackingPath(YoungControls, config);
%transform data
YoungControls = TransformPaths(YoungControls);

%% 1,the degraded LeakyIntegrated model
config.ModelName = "degradedLIFull";
config.NumParams = 5;
Results_degradedLIFull = PerformGroupFit(YoungControls, config);
IC_degradedLIFull = Results_degradedLIFull.IC;


%% 2,the LeakyIntegrate full model
config.ModelName = "LIFull";
config.NumParams = 5;
Results_LIFull = PerformGroupFit(YoungControls, config);
IC_LIFull = Results_LIFull.IC;

%%
ModelNames = {'degradedLIFull', 'LIFull'};

%% Setting colors for using in plots
ColorPattern; 

%% Box Plot of AIC
ICType = "AIC";

newIC_degradedLIFull = getRawIC(IC_degradedLIFull, ICType);
newIC_LIFull = getRawIC(IC_LIFull, ICType);

IC_Cond1 = [newIC_degradedLIFull', newIC_LIFull'];
plotBoxPlot(IC_Cond1, ModelNames, ICType, config);

% Box Plot of BIC 
ICType = "BIC";

newIC_degradedLIFull = getRawIC(IC_degradedLIFull, ICType);
newIC_LIFull = getRawIC(IC_LIFull, ICType);

IC_Cond1 = [newIC_degradedLIFull', newIC_LIFull'];
plotBoxPlot(IC_Cond1, ModelNames, ICType, config);

% Box Plot of NegLogLikelihood
ICType = "NegLogLikelihood";

newIC_degradedLIFull = getRawIC(IC_degradedLIFull, ICType);
newIC_LIFull = getRawIC(IC_LIFull, ICType);

IC_Cond1 = [newIC_degradedLIFull', newIC_LIFull'];
plotBoxPlot(IC_Cond1, ModelNames, ICType, config);

%% function for get the IC data out from the data structure
function icdata=getRawIC(IC, ICType)
    if ICType == "AIC"
        icdata = [];
        for j=1:length(IC)
            icdata = [icdata,IC{j}.aic];
        end
    elseif ICType == "BIC"
        icdata = [];
        for j=1:length(IC)
            icdata = [icdata,IC{j}.bic];
        end
    elseif ICType == "NegLogLikelihood"
        icdata = [];
        for j=1:length(IC)
            icdata = [icdata,IC{j}.negll];
        end      
    else
        error("Choose correct IC type!")
    end
end

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
end
