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

%% BarScatter Plot between HealthyOld and MergedMCI on Fitted Param: Gamma
ParamIndx = 1; 
plotBarScatterOfFittedParam(AllHealthyOldParams, AllMCIParams, ParamIndx, config);

%% BarScatter Plot between HealthyOld and MergedMCI on Fitted Param: g3
ParamIndx = 4; 
plotBarScatterOfFittedParam(AllHealthyOldParams, AllMCIParams, ParamIndx, config);

%% BarScatter Plot between HealthyOld and MergedMCI on Fitted Param: b
ParamIndx = 5; 
plotBarScatterOfFittedParam(AllHealthyOldParams, AllMCIParams, ParamIndx, config);

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

function plotBarScatterOfFittedParam(AllHealthyOldParams, AllMCIParams, ParamIndx, config)
    % AllParams: estimated parameter values
    % ParamIndx
    % config
    
    numConds = 3; %3 is the condition number
    mean_all = zeros(numConds, 2); %2 is the group number
    sem_all = zeros(numConds, 2); %2 is the group number

    ParamName = ["Gamma", "G3", "g2", "g3", "b", "sigma", "nu"];
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
    b(2).FaceColor = colorForMCI; %set color for MCI merged

    hold on
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(mean_all);
    % Get the x coordinate of the bars
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end
    % Plot the errorbars
    errorbar(x',mean_all,sem_all,'k','LineStyle','None', 'LineWidth', 2);
    
    hold on
    % scatter
    for TRIAL_FILTER=1:numConds
        %% extract data
        HealthyOldParam = AllHealthyOldParams{TRIAL_FILTER}(:,ParamIndx);
        MCIParam = AllMCIParams{TRIAL_FILTER}(:,ParamIndx);    

        %plot scatter for HealthyOld
        for i=1:length(HealthyOldParam)
            jitter_value = 0.2*(rand(1)-0.5);
            scatter(x(1,TRIAL_FILTER)+jitter_value, HealthyOldParam(i),30, ...
            "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForHOld, ...
            'MarkerEdgeAlpha', 1, 'MarkerFaceAlpha', 0.3);
        end

        hold on
        %plot scatter for MCI
        for i=1:length(MCIParam)
            jitter_value = 0.2*(rand(1)-0.5);
            scatter(x(2,TRIAL_FILTER)+jitter_value, MCIParam(i),30, ...
            "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colorForMCI, ...
            'MarkerEdgeAlpha', 1, 'MarkerFaceAlpha', 0.3);  
        end

    end    
    
    %Further post-processing the figure
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XLim'        , [0.5, 3.5],...
        'YLim'        , [0, 1.5],...   
        'XTick'       , [1:3],... 
        'XTickLabel'  , {'No Change','No Distal Cue', 'No Optical Flow'},...
        'Ytick'       , [0,0.5,1.0,1.5],...
        'LineWidth'   , .5        );
    ylabel('\gamma');
    legend(b, {'HealthyOld' 'MCIMerged'}, 'Location','northwest', 'NumColumns',2);
    %xlabel('G_1','Interpreter','tex'); ylabel('G_2','Interpreter','tex');
    
    %% save figure
    exportgraphics(f,config.ResultFolder+"/Bar_"+ParamName(ParamIndx)+".png",'Resolution',300);
    exportgraphics(f,config.ResultFolder+"/Bar_"+ParamName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');
end

function [m, s] = getMeanSem(Param)
    m = mean(Param);
    s = std(Param)./sqrt(length(Param));
end