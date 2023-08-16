%% Script that comparing leaky integration with time and distance
% Zilong Ji, UCL, 2022 zilong.ji@ucl.ac.uk
% We calculated different generative angular-linear models with different
% types of source errors to find the candidate for our model (please see
% online methods for details. We chose healthy elderly group to fit our
% models. For each model we calculate the Akaike information criterion
% (AIC) and the Bayesian criterion information (BIC) to select the final
% model taking into account the model complexity.

% Preparing the data
GLAMPI_PrepareBaseConfig

% Preprocessing the data
GLAMPI_PreprocessData

% Instructing the model to run only on healthy elderly
config.ModelSelection = true;

%% leaky integration with time (Model in the paper)
rng('default');

config.ModelName        =   "beta_k_g2_g3_m3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "m3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

LIWithTime_IC  =   HealthyControls.Results.IC;


% leaky integration with distance
config.ModelName        =   "alpha_k_g2_g3_m3_sigma_nu";
config.ParamName        =   ["alpha", "k", "g2", "g3", "m3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

LIWithDist_IC  =   HealthyControls.Results.IC;


%% Preparing output folder
config.ResultFolder = pwd+"/Output/FigS9";

if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Generating color scheme for our paper
ColorPattern; 

% Generating schematic names for our model -> MX.Y where X is total number
% of fitted parameters and Y is variant within a given X
ModelNames = {'LIWithTime', 'LIWithDist'};

%% reformating the IC structure
CondList = ["no change", "no distal cue", "no optical flow", "all"];

% When cond = "all" it will average the values between the conditions
for idx_ = 1:4
    cond = CondList(idx_);

    [LIWithTime_AIC,LIWithTime_BIC, LIWithTime_NLL] = reformatIC(LIWithTime_IC, cond);
    [LIWithDist_AIC,LIWithDist_BIC, LIWithDist_NLL] = reformatIC(LIWithDist_IC, cond);

    % Box Plot for AICs
    ICType = "AIC";
    All_AIC = [LIWithTime_AIC, LIWithDist_AIC];

    plotBoxPlot(All_AIC, ModelNames, ICType, config, cond);
    plotErrorPlot(All_AIC, ModelNames, ICType, config, cond);
    
    % Box Plot for BIC
    ICType = "BIC";
    All_BIC = [LIWithTime_BIC, LIWithDist_BIC];

    plotBoxPlot(All_BIC, ModelNames, ICType, config, cond);
    plotErrorPlot(All_BIC, ModelNames, ICType, config, cond);
    
    % Box Plot for NLL
    ICType = "NegLogLikelihood";
    All_NLL = [LIWithTime_NLL, LIWithDist_NLL];

    plotBoxPlot(All_NLL, ModelNames, ICType, config, cond);
    plotErrorPlot(All_NLL, ModelNames, ICType, config, cond);

end

% Cleaning dummy variable to allow other scripts to run on all groups
if (isfield(config, 'ModelSelection') == 1)
    config = rmfield(config,"ModelSelection");
end

function plotBoxPlot(data, ModelNames, ICType, config, cond)

    f = figure('visible','off','Position', [100 100 600 400]);
    
    % Parameters set for controlling visual output
    num_boxplots = size(data,2);
    box_lineWidth = 0.5;
    whisker_value = 1.5;
    box_widths_value = 0.6;
    box_color_transparency = 0.5; 
    median_lineWidth = 2;
    median_color = 'k';
    scatter_jitter_value = 0.1;
    scatter_markerSize=10;
    scatter_marker_edgeColor = 'k';
    scatter_marker_edgeWidth = 0.5;
    scatter_color_transparency = 0.7; 

    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)
    
    N = size(config.color_scheme_npg,1)-2;
    color_scheme = config.color_scheme_npg(1:N,:);
    
    bp = boxplot(data, 'whisker',whisker_value,'symbol','', ... 
                'color','k', ...
                'Notch','on',...
                'widths',box_widths_value,...
                'labels', ModelNames);
    
    set(bp,'linewidth',box_lineWidth);

    if ICType=="AIC"
        ylabel('Akaike information criterion (AIC)');
    elseif ICType=="BIC"
        ylabel('Bayesian information Criterion (BIC)');
    elseif ICType=="NegLogLikelihood"
        ylabel('Negative Loglikelihood');
    else
        error("Choose correct IC type!");
    end
    
    %% Boxplot visual changes
    h = findobj(gca,'Tag','Box');
    for i = length(h):-1:1
        c_idx = mod(i,N)+1;
        patch(get(h(i),'XData'),get(h(i),'YData'),color_scheme(c_idx,:),'FaceAlpha',box_color_transparency);
    end
    % Sending patch to back of the figure so that median can be drawn on top of it
    set(gca,'children',flipud(get(gca,'children'))) 

    %% Median visual changes
    h=findobj(gca,'tag','Median');
    for i = 1:length(h)
        h(i).LineWidth = median_lineWidth;
        h(i).Color = median_color;
    end

    %% Scatter plot for data and mean
    num_points = size(data,1);

    for i=1:size(data,2)
        hold on
        c_idx = mod(num_boxplots-i+1,N)+1;

        x = i*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, data(:,i), scatter_markerSize, ...
                'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor', color_scheme(c_idx,:),...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 
        hold on
        scatter(i, mean(data(:,i)), 4*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);
    end

    %% Figure post-processing  
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , 1:1:100,... 
        'XLim'        , [0.5,length(ModelNames)+0.5],...
        'LineWidth'   , .5        );
    
    %% Export figure
    exportgraphics(f,config.ResultFolder+"/"+cond+"_"+ICType+".png",'Resolution',300);
    exportgraphics(f,config.ResultFolder+"/"+cond+"_"+ICType+".pdf",'Resolution',300,'ContentType','vector');
end

function plotErrorPlot(data, ModelNames, ICType, config, cond)
    
    f = figure('visible','off','Position', [100 100 400 400]);

    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)
   
    mean_ = mean(data,1, 'omitnan');
    std_ = std(data,0, 1, 'omitnan');
    sem_ = std_./sqrt(size(data,1));
    
    num_boxplots = size(data,2);

    bar(1:1:num_boxplots, mean_, 'FaceColor', config.color_scheme_npg(2,:), 'FaceAlpha', 0.5);

    hold on

    h=errorbar(1:1:num_boxplots, mean_, sem_, 'LineStyle','none', 'Color', 'k','linewidth', 3);
    h.CapSize = 20;

    hold on
    scatter(1:1:num_boxplots, mean_, 100, 'd',...
        'filled','MarkerEdgeColor','k', ...
        'MarkerFaceColor','r', ...
        'LineWidth',0.5);

    if ICType=="AIC"
        ylabel('Akaike information criterion (AIC)');
    elseif ICType=="BIC"
        ylabel('Bayesian Inference Criterion (BIC)');
    elseif ICType=="NegLogLikelihood"
        ylabel('Negative Loglikelihood');
    else
        error("Choose correct IC type!");
    end

    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , 1:1:num_boxplots,... 
        'XTickLabel'  , ModelNames,...
        'XLim'        , [0.5,num_boxplots+0.5],...
        'LineWidth'   , .5        );

    %% Export figures
    exportgraphics(f,config.ResultFolder+"/ErrorBar_"+cond+"_"+ICType+".png",'Resolution',300);
    exportgraphics(f,config.ResultFolder+"/ErrorBar_"+cond+"_"+ICType+".pdf",'Resolution',300,'ContentType','vector');
end