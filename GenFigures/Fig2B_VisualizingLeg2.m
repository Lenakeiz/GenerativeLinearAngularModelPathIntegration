%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default'); %for code reproducibility

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
YoungControls = TransformPaths(YoungControls);
savefolder = pwd + "/Output/";

%% setting the configuration
config.UseGlobalSearch = true;
resultfolder = savefolder+"PaperFigs/Fig2B";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% Model fitting using our base model
%no consider of the rotation gain factor in leg2, i.e., gamma, G3=1, g2=1, g3, k3, sigma, nu. #params=6
config.ModelName = "BaseModel";
config.NumParams = 5;
[AllYoungParams, AllYoungX, AllYoungDX, AllYoungTheta, AllYoungIC] = getResultsAllConditions(YoungControls, config);

%% Setting colors for using in plots
ColorPattern; 

%% Plot the leg 1 distribution
plotLeg2(AllYoungParams, AllYoungDX, config)

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

function plotLeg2(AllParams, AllDX, config)
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
        background_color = config.color_scheme_npg(:,:);
        forground_color = [0.8500 0.3250 0.0980];
        numcolors = size(background_color,1);

        %% extract data
        ParamValues = AllParams{TRIAL_FILTER};
        DX = AllDX{TRIAL_FILTER};
        subjectSize = size(DX,2);

        all_Leg2 = []; all_Leg2_prime = [];
        
        for subj=1:subjectSize %for each subject
            paramX = ParamValues(subj,:);
            subjDX = DX{subj};
            %extract gamma and G3
            gamma = paramX(1); G3 = paramX(2);
            
            sampleSize = size(subjDX,2);

            Leg2 = zeros(1,sampleSize);
            Leg2_prime = zeros(1,sampleSize);
            
            %for each subject choose one color by looping over color_scheme_npg
            color_idx = mod(subj,numcolors)+1;
            cm = background_color(color_idx,:);
            %add jitter to making better plot
            jitter_value = 0.6*(rand(1)-0.5);
            for tr_id = 1:sampleSize %for each trial
                leg2 = subjDX{tr_id}(1); Leg2(tr_id)=leg2;
                leg2_prime = gamma*G3*leg2; Leg2_prime(tr_id)=leg2_prime;
                pt = plot([1+jitter_value, 2+jitter_value], [leg2, leg2_prime], '-', 'Color', cm ,'linewidth',2);            
                %transparancy
                pt.Color = [pt.Color 0.05];
                hold on
                %scatter
                st = scatter([1+jitter_value, 2+jitter_value], [leg2, leg2_prime], 15, ...
                    "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cm, ...
                    'MarkerEdgeAlpha', 0.6, 'MarkerFaceAlpha', 0.3);
                hold on
            end
            all_Leg2 = [all_Leg2,Leg2];
            all_Leg2_prime = [all_Leg2_prime,Leg2_prime];
        end

        all_Leg1_mean = mean(all_Leg2);
        all_Leg1_prime_mean = mean(all_Leg2_prime);
        pt = plot([1 2], [all_Leg1_mean all_Leg1_prime_mean], '-d', 'Color', forground_color, 'linewidth',8, 'MarkerSize', 8);
        pt.Color = [pt.Color 0.8];

        ylabel('Distance (in meters)');
        %Further post-processing the figure
        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.1 .1 .1], ...
            'YColor'      , [.1 .1 .1], ...
            'XTick'       , [1,2],... 
            'XTickLabel'  , {'Ideal distance','Mental distance'},...
            'YLim'        , [min([all_Leg2,all_Leg2_prime])-0.1,max([all_Leg2,all_Leg2_prime])+0.1],...
            'XLim'        , [0.5, 2.5],...
            'LineWidth'   , .5        );
            %'Ytick'       , [1:ceil(max([all_Leg2,all_Leg2_prime]))],...
            %'YLim'        , [0,ceil(max([all_Leg2,all_Leg2_prime]))],...
        title('Outbound leg 2');
        %% save figure
        exportgraphics(f,config.ResultFolder+"/Leg2"+conditionName{TRIAL_FILTER}+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/Leg2"+conditionName{TRIAL_FILTER}+".pdf",'Resolution',300, 'ContentType','vector');
    end
end



