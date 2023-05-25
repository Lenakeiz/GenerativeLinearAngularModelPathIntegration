%% Script that fit different types of the generative linear angular model considering differnt types of error - output is in Fig. 2, Fig S2, Fig S4
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

%% 1, sigma nu Model --> egocentric noise only model
config.ModelName        =   "sigma_nu";
config.ParamName        =   ["sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

sigma_nu_IC             =   HealthyControls.Results.IC;

%% 2.1, beta, k, sigma, nu Model --> only encoding errors in distance
config.ModelName        =   "beta_k_sigma_nu";
config.ParamName        =   ["beta", "k", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

beta_k_sigma_nu_IC        =   HealthyControls.Results.IC;

%% 2.2, g2, sigma, nu Model --> only encoding errors in distance
config.ModelName        =   "g2_sigma_nu";
config.ParamName        =   ["g2", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

g2_sigma_nu_IC        =   HealthyControls.Results.IC;

%% 2.3, g3, sigma, nu Model --> only encoding errors in distance
config.ModelName        =   "g3_sigma_nu";
config.ParamName        =   ["g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

g3_sigma_nu_IC        =   HealthyControls.Results.IC;

%% 2.4, m3, sigma, nu Model --> only encoding errors in distance
config.ModelName        =   "m3_sigma_nu";
config.ParamName        =   ["m3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

m3_sigma_nu_IC        =   HealthyControls.Results.IC;

%% 3.1 beta,k,g2,sigma,nu model --> encoding errors in both distance and angle
config.ModelName        =   "beta_k_g2_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

beta_k_g2_sigma_nu_IC        =   HealthyControls.Results.IC;

%% 3.2 beta,k,g3,sigma,nu model --> encoding errors in both distance and angle
config.ModelName        =   "beta_k_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

beta_k_g3_sigma_nu_IC        =   HealthyControls.Results.IC;

%% 3.3 beta,k,m3,sigma,nu model --> encoding errors in both distance and angle
config.ModelName        =   "beta_k_m3_sigma_nu";
config.ParamName        =   ["beta", "k", "m3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

beta_k_m3_sigma_nu_IC        =   HealthyControls.Results.IC;


%% 3.4, g2, g3, sigma, nu Model --> encoding errors and production errors in angle 
config.ModelName        =   "g2_g3_sigma_nu";
config.ParamName        =   ["g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

g2_g3_sigma_nu_IC       =   HealthyControls.Results.IC;

%% 3.5, g2, m3, sigma, nu Model --> encoding errors and production errors in angle 
config.ModelName        =   "g2_m3_sigma_nu";
config.ParamName        =   ["g2", "m3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

g2_m3_sigma_nu_IC       =   HealthyControls.Results.IC;

%% 3.6, g3, m3, sigma, nu Model --> encoding errors and production errors in angle 
config.ModelName        =   "g3_m3_sigma_nu";
config.ParamName        =   ["g3", "m3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

g3_m3_sigma_nu_IC       =   HealthyControls.Results.IC;

%% 4.1 beta, k, g2, g3, sigma, nu Model --> encoding errors in distance and angles and production error in angle
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

beta_k_g2_g3_sigma_nu_IC  =   HealthyControls.Results.IC;

%% 4.2 beta, k, g2, m3, sigma, nu Model --> encoding errors in distance and angles and production error in angle
config.ModelName        =   "beta_k_g2_m3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "m3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

beta_k_g2_m3_sigma_nu_IC  =   HealthyControls.Results.IC;

%% 4.3 beta, k, g3, m3, sigma, nu Model --> encoding errors in distance and angles and production error in angle
config.ModelName        =   "beta_k_g3_m3_sigma_nu";
config.ParamName        =   ["beta", "k", "g3", "m3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

beta_k_g3_m3_sigma_nu_IC  =   HealthyControls.Results.IC;

%% 4.4 g2, g3, m3, sigma, nu Model --> encoding errors in distance and angles and production error in angle
config.ModelName        =   "g2_g3_m3_sigma_nu";
config.ParamName        =   ["g2", "g3", "m3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

g2_g3_m3_sigma_nu_IC  =   HealthyControls.Results.IC;

%% 5.1 beta, k, g2, g3, m3, sigma, nu Model --> encoding errors in distance and angles and production error in angle
config.ModelName        =   "beta_k_g2_g3_m3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "m3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

beta_k_g2_g3_m3_sigma_nu_IC  =   HealthyControls.Results.IC;


%% beta1 beta2 Model --> vector addition model with separate gains for the outbound path segments, follows Harootonian et al., 2020 paper
config.ModelName        =   "beta1_beta2_g2_g3_sigma_nu";
config.ParamName        =   ["beta1", "beta2", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI

beta1_beta2_g2_g3_sigma_nu_IC  =   HealthyControls.Results.IC;

% Preparing output folder
config.ResultFolder = pwd+"/Output/ModelSelection/GenerativeLinearAngularModel";

if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Generating color scheme for our paper
ColorPattern; 

% Generating schematic names for our model -> MX.Y where X is total number
% of fitted parameters and Y is variant within a given X
ModelNames = {'M2', ...
              'M3.1', 'M3.2', 'M3.3', 'M3.4',...
              'M4.1', 'M4.2','M4.3','M4.4','M4.5','M4.6',...
              'M5.1', 'M5.2','M5.3','M5.4',...
              'M6'};

%% reformating the IC structure
CondList = ["no change", "no distal cue", "no optical flow", "all"];

% When cond = "all" it will average the values between the conditions
for idx_ = 1:4
    cond = CondList(idx_);

    [sigma_nu_AIC,sigma_nu_BIC, sigma_nu_NLL] = reformatIC(sigma_nu_IC, cond);
    [beta_k_sigma_nu_AIC,beta_k_sigma_nu_BIC,beta_k_sigma_nu_NLL] = reformatIC(beta_k_sigma_nu_IC, cond);
    [g2_sigma_nu_AIC,g2_sigma_nu_BIC, g2_sigma_nu_NLL] = reformatIC(g2_sigma_nu_IC, cond);
    [g3_sigma_nu_AIC,g3_sigma_nu_BIC, g3_sigma_nu_NLL] = reformatIC(g3_sigma_nu_IC, cond);
    [m3_sigma_nu_AIC,m3_sigma_nu_BIC, m3_sigma_nu_NLL] = reformatIC(m3_sigma_nu_IC, cond);
    [beta_k_g2_sigma_nu_AIC,beta_k_g2_sigma_nu_BIC, beta_k_g2_sigma_nu_NLL] = reformatIC(beta_k_g2_sigma_nu_IC, cond);
    [beta_k_g3_sigma_nu_AIC,beta_k_g3_sigma_nu_BIC, beta_k_g3_sigma_nu_NLL] = reformatIC(beta_k_g3_sigma_nu_IC, cond);
    [beta_k_m3_sigma_nu_AIC,beta_k_m3_sigma_nu_BIC, beta_k_m3_sigma_nu_NLL] = reformatIC(beta_k_m3_sigma_nu_IC, cond);
    [g2_g3_sigma_nu_AIC,g2_g3_sigma_nu_BIC, g2_g3_sigma_nu_NLL] = reformatIC(g2_g3_sigma_nu_IC, cond);
    [g2_m3_sigma_nu_AIC,g2_m3_sigma_nu_BIC, g2_m3_sigma_nu_NLL] = reformatIC(g2_m3_sigma_nu_IC, cond);
    [g3_m3_sigma_nu_AIC,g3_m3_sigma_nu_BIC, g3_m3_sigma_nu_NLL] = reformatIC(g3_m3_sigma_nu_IC, cond);
    [beta_k_g2_g3_sigma_nu_AIC,beta_k_g2_g3_sigma_nu_BIC, beta_k_g2_g3_sigma_nu_NLL] = reformatIC(beta_k_g2_g3_sigma_nu_IC, cond);
    [beta_k_g2_m3_sigma_nu_AIC,beta_k_g2_m3_sigma_nu_BIC, beta_k_g2_m3_sigma_nu_NLL] = reformatIC(beta_k_g2_m3_sigma_nu_IC, cond);
    [beta_k_g3_m3_sigma_nu_AIC,beta_k_g3_m3_sigma_nu_BIC, beta_k_g3_m3_sigma_nu_NLL] = reformatIC(beta_k_g3_m3_sigma_nu_IC, cond);
    [g2_g3_m3_sigma_nu_AIC,g2_g3_m3_sigma_nu_BIC, g2_g3_m3_sigma_nu_NLL] = reformatIC(g2_g3_m3_sigma_nu_IC, cond);
    [beta_k_g2_g3_m3_sigma_nu_AIC,beta_k_g2_g3_m3_sigma_nu_BIC, beta_k_g2_g3_m3_sigma_nu_NLL] = reformatIC(beta_k_g2_g3_m3_sigma_nu_IC, cond);
    
    % Box Plot for AICs
    ICType = "AIC";
    All_AIC = [sigma_nu_AIC,...              %2
               beta_k_sigma_nu_AIC,...       %3.1
               g2_sigma_nu_AIC,...           %3.2
               g3_sigma_nu_AIC,...           %3.3
               m3_sigma_nu_AIC,...           %3.4
               beta_k_g2_sigma_nu_AIC,...    %4.1
               beta_k_g3_sigma_nu_AIC,...    %4.2
               beta_k_m3_sigma_nu_AIC,...    %4.3
               g2_g3_sigma_nu_AIC,...        %4.4
               g2_m3_sigma_nu_AIC,...        %4.5
               g3_m3_sigma_nu_AIC,...        %4.6
               beta_k_g2_g3_sigma_nu_AIC,... %5.1
               beta_k_g2_m3_sigma_nu_AIC,... %5.2
               beta_k_g3_m3_sigma_nu_AIC,... %5.3
               g2_g3_m3_sigma_nu_AIC,...     %5.4
               beta_k_g2_g3_m3_sigma_nu_AIC]; %6.1

    plotBoxPlot(All_AIC, ModelNames, ICType, config, cond);
    plotErrorPlot(All_AIC, ModelNames, ICType, config, cond);
    
    % Box Plot for BIC
    ICType = "BIC";
    All_BIC = [sigma_nu_BIC,...              %2
               beta_k_sigma_nu_BIC,...       %3.1
               g2_sigma_nu_BIC,...           %3.2
               g3_sigma_nu_BIC,...           %3.3
               m3_sigma_nu_BIC,...           %3.4
               beta_k_g2_sigma_nu_BIC,...    %4.1
               beta_k_g3_sigma_nu_BIC,...    %4.2
               beta_k_m3_sigma_nu_BIC,...    %4.3
               g2_g3_sigma_nu_BIC,...        %4.4
               g2_m3_sigma_nu_BIC,...        %4.5
               g3_m3_sigma_nu_BIC,...        %4.6
               beta_k_g2_g3_sigma_nu_BIC,... %5.1
               beta_k_g2_m3_sigma_nu_BIC,... %5.2
               beta_k_g3_m3_sigma_nu_BIC,... %5.3
               g2_g3_m3_sigma_nu_BIC,...     %5.4
               beta_k_g2_g3_m3_sigma_nu_BIC];%6.1

    plotBoxPlot(All_BIC, ModelNames, ICType, config, cond);
    plotErrorPlot(All_BIC, ModelNames, ICType, config, cond);
    
    % Box Plot for NLL
    ICType = "NegLogLikelihood";
    All_NLL = [sigma_nu_NLL,...              %1
               beta_k_sigma_nu_NLL,...       %2.1
               g2_sigma_nu_NLL,...           %2.2
               g3_sigma_nu_NLL,...           %2.3
               m3_sigma_nu_NLL,...           %2.4
               beta_k_g2_sigma_nu_NLL,...    %3.1
               beta_k_g3_sigma_nu_NLL,...    %3.2
               beta_k_m3_sigma_nu_NLL,...    %3.3
               g2_g3_sigma_nu_NLL,...        %3.4
               g2_m3_sigma_nu_NLL,...        %3.5
               g3_m3_sigma_nu_NLL,...        %3.6
               beta_k_g2_g3_sigma_nu_NLL,... %4.1
               beta_k_g2_m3_sigma_nu_NLL,... %4.2
               beta_k_g3_m3_sigma_nu_NLL,... %4.3
               g2_g3_m3_sigma_nu_NLL,...     %4.4
               beta_k_g2_g3_m3_sigma_nu_NLL];%5.1

    plotBoxPlot(All_NLL, ModelNames, ICType, config, cond);
    plotErrorPlot(All_NLL, ModelNames, ICType, config, cond);

end

% Final cleanup to leave workspace as the end of the Preprocessing stage.
% Remove if you want to take a look at the output data.
clearvars -except config YoungControls HealthyControls MCINeg MCIPos MCIUnk
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
    
    f = figure('visible','off','Position', [100 100 1200 300]);

    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)
   
    mean_ = mean(data,1, 'omitnan');
    std_ = std(data,0, 1, 'omitnan');
    sem_ = std_./sqrt(size(data,1));
    
    num_boxplots = size(data,2);
    h=errorbar(1:1:num_boxplots, mean_, sem_, 'LineStyle','none', 'Color', 'k','linewidth', 3);
    h.CapSize = 20;

    hold on

    scatter(1:1:num_boxplots, mean_, 50, 'd',...
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