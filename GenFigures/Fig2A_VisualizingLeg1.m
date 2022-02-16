%% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
YoungControls = TransformPaths(YoungControls);
savefolder = pwd + "/Output/";

%% setting the configuration
config.UseGlobalSearch = true;
resultfolder = savefolder+"PaperFigs/Fig2A";
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

%% Plot the leg 1 distribution
plotLeg1(AllYoungParams, AllYoungDX, resultfolder)
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

function plotLeg1(AllParams, AllDX, resultfolder)
    %%plot leg 1 changes
    % AllParams: estimated parameter values
    % AllDX: is a cell structure containing the segment of each trial
    
    numConds = length(AllParams); %3
    conditionName = {'No Change', 'No Distal Cue', 'No Optical Flow'};
    for TRIAL_FILTER=1:numConds
        
        %% set figure info
        f = figure('visible','off','Position', [100 100 600 400]);
        %%% Font type and size setting %%%
        % Using Arial as default because all journals normally require the font to
        % be either Arial or Helvetica
        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)     

        %%% Color definition %%%
        color_scheme_npg = [0         0.6275    0.5294; ... 
                            0.9020    0.2941    0.2078; ...
                            0.3020    0.7333    0.8353; ...
                            0.9529    0.6078    0.4980; ...                     
                            0.2353    0.3294    0.5333; ...
                            0.5176    0.5686    0.7059; ...
                            0.5686    0.8196    0.7608; ...
                            0.8627         0         0; ...
                            0.4941    0.3804    0.2824; ...
                            0.6902    0.6118    0.5216 ];   
        numcolors = size(color_scheme_npg,1);

        %% extract data
        ParamValues = AllParams{TRIAL_FILTER};
        DX = AllDX{TRIAL_FILTER};
        subjectSize = size(DX,2);

        all_Leg1 = []; all_Leg1_prime = [];
        
        for subj=1:subjectSize %for each subject
            paramX = ParamValues(subj,:);
            subjDX = DX{subj};
            %extract gamma and G3
            gamma = paramX(1); G3 = paramX(2);
            
            sampleSize = size(subjTHETADX,2);

            Leg1 = zeros(1,sampleSize);
            Leg1_prime = zeros(1,sampleSize);
            
            %for each subject randomly choose one color from color_scheme_npg
            color_idx = randsample(numcolors,1);
            cm = color_scheme_npg(color_idx,:);
            for tr_id = 1:sampleSize %for each trial
                leg1 = subjDX{tr_id}(1); Leg1(tr_id)=leg1;
                leg1_prime = gamma^2*G3*leg1; Leg1_prime(tr_id)=leg1_prime;
                pt = plot([1, 2], [leg1, leg1_prime], '-o', "Color", cm, 'linewidth',2, 'MarkerSize', 0.5);            
                %transparancy
                pt.Color = [pt.Color 0.05];
                hold on
            end
            all_Leg1 = [all_Leg1,Leg1];
            all_Leg1_prime = [all_Leg1_prime,Leg1_prime];
        end

        all_Leg1_mean = mean(all_Leg1);
        all_Leg1_prime_mean = mean(all_Leg1_prime);
        plot([1 2], [all_Leg1_mean all_Leg1_prime_mean], '-dk', 'linewidth',5, 'MarkerSize', 5)

        ylabel('Distance (meters)');
        %Further post-processing the figure
        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.1 .1 .1], ...
            'YColor'      , [.1 .1 .1], ...
            'XTick'       , [1,2],... 
            'XTickLabel'  , {'Physical','Mental'},...
            'XLim'        , [0.5, 2.5],...
            'YLim'        , [0,ceil(max(all_Leg1))],...
            'LineWidth'   , .5        );
        title(conditionName{TRIAL_FILTER});


    end
    %% save figure
    exportgraphics(f,resultfolder+"/Leg1Change.png",'Resolution',300);
end



