%% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
YoungControls = TransformPaths(YoungControls);
savefolder = pwd + "/Output/";

%% setting the configuration
config.UseGlobalSearch = true;
resultfolder = savefolder+"PaperFigs/Fig3A";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% Model fitting using our base model
%no consider of the rotation gain factor in leg2, i.e., gamma, G3=1, g2=1, g3, k3, sigma, nu. #params=6
config.ModelName = "G1G2";
config.NumParams = 6;
[AllYoungParams, AllYoungX, AllYoungDX, AllYoungTheta, AllYoungIC] = getResultsAllConditions(YoungControls, config);

%% Setting colors for using in plots
ColorPattern; 

%% Plot the leg 1 distribution
plotRegressionBetweenG1G2Square(AllYoungParams, config)

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

function plotRegressionBetweenG1G2Square(AllParams, config)
    %%plot leg 1 changes
    % AllParams: estimated parameter values
    % AllDX: is a cell structure containing the segment of each trial
    
    numConds = length(AllParams); %3
    conditionName = {"No Change", "No Distal Cue", "No Optical Flow"};
    for TRIAL_FILTER=1:numConds
        
        %% set figure info
        f = figure('visible','off','Position', [100 100 700 500]);

        %%% Font type and size setting %%%
        % Using Arial as default because all journals normally require the font to
        % be either Arial or Helvetica
        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)     

        %%% Color definition %%%
        color_scheme_npg = config.color_scheme_npg;

        %% extract data
        ParamValues = AllParams{TRIAL_FILTER};
        subjectSize = size(ParamValues,1);
        all_G1 = zeros(1,subjectSize); all_G2 = zeros(1,subjectSize);
        for subj=1:subjectSize %for each subject
            paramX = ParamValues(subj,:);
            %extract G1 and G2
            G1 = paramX(1); all_G1(subj)=G1;
            G2 = paramX(2); all_G2(subj)=G2;
        end

        md1 = fitlm(all_G1,all_G2.^2, 'linear'); %returns a linear regression between G1 and G2^2
        h = plot(md1, 'LineWidth', 3);
        set(h(1), 'Color', 'r', 'Marker', '.', 'MarkerSize', 20); %set data point
        set(h(2), 'Color', 'k', 'LineWidth', 3); %set the regression line
        set(h(3), 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.'); %set the upper confidence line
        set(h(4), 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.'); %set the lower confidence line

        %Further post-processing the figure
        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.1 .1 .1], ...
            'YColor'      , [.1 .1 .1], ...
            'XLim'        , [0, 1.5],...
            'YLim'        , [0, 1.5],...        
            'LineWidth'   , .5        );
            %'XTick'       , [0,0.5,1.0],... 
            %'Ytick'       , [0,0.5,1.0],...
        xlabel('G_1','Interpreter','tex'); ylabel('G_2^2','Interpreter','tex');
        
        %title(conditionName{TRIAL_FILTER});
        [R,P]=corrcoef(all_G1,all_G2.^2);
        title(conditionName{TRIAL_FILTER}+" Person Corrcoef="+num2str(round(R(1,2),4))+" P="+num2str(round(P(1,2),4)));
        %% save figure
        exportgraphics(f,config.ResultFolder+"/G1_G2Square_"+conditionName{TRIAL_FILTER}+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/G1_G2Square_"+conditionName{TRIAL_FILTER}+".pdf",'Resolution',300, 'ContentType','vector');
    end
end



