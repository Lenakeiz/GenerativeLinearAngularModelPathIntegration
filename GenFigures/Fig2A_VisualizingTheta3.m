%% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
YoungControls = TransformPaths(YoungControls);
savefolder = pwd + "/Output/";

%% setting the configuration
config.UseGlobalSearch = true;
resultfolder = savefolder+"PaperFigs/Fig2C";
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
plotTheta3(AllYoungParams, AllYoungDX, AllYoungTheta, config)

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

function plotTheta3(AllParams, AllDX, AllTheta, config)
    %%plot leg 1 changes
    % AllParams: estimated parameter values
    % AllDX: is a cell structure containing the segment of each trial
    
    numConds = length(AllParams); %3
    conditionName = {"No Change", "No Distal Cue", "No Optical Flow"};
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
        background_color = config.color_scheme_npg(:,:);
        forground_color = [0.8500 0.3250 0.0980];
        numcolors = size(background_color,1);

        %% extract data
        ParamValues = AllParams{TRIAL_FILTER};
        DX = AllDX{TRIAL_FILTER};
        THETADX = AllTheta{TRIAL_FILTER};
        subjectSize = size(DX,2);

        all_alpha = []; all_theta3_prime = [];
        
        for subj=1:subjectSize %for each subject
            paramX = ParamValues(subj,:);
            subjDX = DX{subj};
            subjTHETADX = THETADX{subj};
            %extract parameters
            [gamma,G3,g2,g3,b,sigma,nu] = deal(paramX(1),paramX(2),paramX(3),paramX(4),paramX(5),paramX(6),paramX(7));
            
            sampleSize = size(subjDX,2);

            Alpha = zeros(1,sampleSize);
            Theta3_prime = zeros(1,sampleSize);
            
            %for each subject choose one color by looping over color_scheme_npg
            color_idx = mod(subj,numcolors)+1;
            cm = background_color(color_idx,:);
            %add jitter to making better plot
            jitter_value = 0.3*(rand(1)-0.5);

            for tr_id = 1:sampleSize %for each trial
                theta3 = subjTHETADX{tr_id}(3); 
                %Theta3(tr_id) = theta3;
                
                %leg length
                l1 = subjDX{tr_id}(1); l2 = subjDX{tr_id}(2);
                %angle value
                theta2 = subjTHETADX{tr_id}(2); theta2_prime = g2*theta2;

                %mental leg 1 end point:
                G1 = gamma^2*G3; men_p1 = [G1*l1, 0];
                
                %mental leg 2 end point:
                G2 = gamma*G3; 
                men_p2 = [G1*l1+G2*l2*cos(theta2_prime), G2*l2*sin(theta2_prime)];

                %for mental leg3 end point:
                %calculate turn angle of accurate mental vector 3
                vec1 = men_p2-men_p1; vec2 = [0,0]-men_p2;
                alpha = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
                %transfer from degree to radians
                alpha = deg2rad(alpha);
                %wrap to (0,2pi)
                alpha = mod(alpha, 2*pi);
                Alpha(tr_id)=alpha;                

                %considering execution turn error
                theta3_prime = g3*alpha+b;
                %also wrap to (0,2pi)
                theta3_prime = mod(theta3_prime, 2*pi);     
                Theta3_prime(tr_id) = theta3_prime;

                pt = plot([1+jitter_value, 2+jitter_value], [alpha, theta3], '-', 'Color', cm ,'linewidth',2);            
                %transparancy
                pt.Color = [pt.Color 0.05];
                hold on
                %scatter
                st = scatter([1+jitter_value, 2+jitter_value], [alpha, theta3], 15, ...
                    "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cm, ...
                    'MarkerEdgeAlpha', 0.6, 'MarkerFaceAlpha', 0.3);
                hold on
            end

            all_alpha = [all_alpha,Alpha];
            all_theta3_prime = [all_theta3_prime,Theta3_prime];  
        end

        all_alpha_mean = mean(all_alpha);
        all_theta3_prime_mean = mean(all_theta3_prime);
        pt = plot([1 2], [all_alpha_mean all_theta3_prime_mean], '-d', 'Color', forground_color, 'linewidth',8, 'MarkerSize', 8);
        pt.Color = [pt.Color 0.8];

        ylabel('Distance (meters)');
        %Further post-processing the figure
        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.1 .1 .1], ...
            'YColor'      , [.1 .1 .1], ...
            'XTick'       , [1,2],... 
            'XTickLabel'  , {'Mental encoded','Mental executed'},...
            'XLim'        , [0.6, 2.4],...
            'YLim'        , [0,ceil(max([all_alpha,all_theta3_prime]))],...
            'LineWidth'   , .5        );
        title(conditionName{TRIAL_FILTER});
        %% save figure
        exportgraphics(f,config.ResultFolder+"/Theta3Change"+conditionName{TRIAL_FILTER}+".png",'Resolution',300);
    end
end



