%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default');

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
HealthyControls = TransformPaths(HealthyControls);
MCIPos          = TransformPaths(MCIPos);
MCINeg          = TransformPaths(MCINeg);
Unknown         = TransformPaths(Unknown);
savefolder = pwd + "/Output/";

%% setting the configuration
resultfolder = savefolder+"PaperFigs/Fig5B";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

ColorPattern;

clear savefolder resultfolder
%% Aggregating the MCI sample
AggregateMCI;

%% Calculating our model, but returning only the informations about the trials
%no consider of the rotation gain factor in leg2, i.e., gamma, G3=1, g2=1, g3, k3, sigma, nu. #params=6
%"VariableNames", {'Gamma','g','b','sigma','nu'});
config.UseGlobalSearch = true;
config.ModelName = "BaseModel";
config.NumParams = 5;

% Model fitting for HealthyOld data
[~, AllHealthyOldX, AllHealthyOldDX, ~, ~] = getResultsAllConditions(HealthyControls, config);
[~, AllMCIAgg, AllMCIAggDX, ~, ~]          = getResultsAllConditions(MCIAll, config);
[~, AllMCINegX, AllMCINegDX, ~, ~]         = getResultsAllConditions(MCINeg, config);
[~, AllMCIPosX, AllMCIPosDX, ~, ~]         = getResultsAllConditions(MCIPos, config);

%% Calculating l3 - l3 hat (walked length - real return length)

[L3_L3Hat_HC]     = getReturnAndActualDistance(AllHealthyOldX, AllHealthyOldDX);
[L3_L3Hat_MCIAll] = getReturnAndActualDistance(AllMCIAgg, AllMCIAggDX);
[L3_L3Hat_MCINeg] = getReturnAndActualDistance(AllMCINegX, AllMCINegDX);
[L3_L3Hat_MCIPos] = getReturnAndActualDistance(AllMCIPosX, AllMCIPosDX);

%% Creating the roc curves for MCI all vs HC
close all;

% set figure info
%f = figure('visible','off','Position', [100 100 1000 500]);
f = figure('visible','off','Position', [100 100 500 500]);
%%% Font type and size setting %%%
% Using Arial as default because all journals normally require the font to
% be either Arial or Helvetica
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12)

hold on;

AUC = plotROCCurve(L3_L3Hat_HC, L3_L3Hat_MCIAll, {'L3_L3Hat'}, 'L3_L3Hat_MCIvsHC','HC', 'MCI', config.color_scheme_npg(4,:), config);
legendText = ['AUC(l3-$\hat{l3}$) = ' , num2str(round(AUC.Value(1),2),2)];

hold off;
    
ll = legend('Location','southeast','Interpreter','latex');
ll.String = legendText;
ll.FontSize = 12;

ylabel('True Positive Rate');
xlabel('False Positive Rate')
title('MCI Merged / Healthy Old controls');

%Further post-processing the figure
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...
    'LineWidth'   , .5        );

axis square;

exportgraphics(f,config.ResultFolder+"/ROCL3_L3Hat_MCIvsHC"+".png",'Resolution',300);
exportgraphics(f,config.ResultFolder+"/ROCL3_L3Hat_MCIvsHC"+".pdf",'Resolution',300, 'ContentType','vector');

clear parametersName filesName i colors f ll legendText 

%% Creating the roc curves for MCI all vs HC
close all;



% set figure info
%f = figure('visible','off','Position', [100 100 1000 500]);
f = figure('visible','off','Position', [100 100 500 500]);
%%% Font type and size setting %%%
% Using Arial as default because all journals normally require the font to
% be either Arial or Helvetica
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12)

hold on;

AUC = plotROCCurve(L3_L3Hat_MCINeg, L3_L3Hat_MCIPos, {'L3_L3Hat'}, 'L3_L3Hat_MCIposvsHCneg','MCIneg', 'MCIpos', config.color_scheme_npg(4,:), config);
legendText = ['AUC(l3-$\hat{l3}$) = ' , num2str(round(AUC.Value(1),2),2)];

hold off;
    
ll = legend('Location','southeast','Interpreter','latex');
ll.String = legendText;
ll.FontSize = 12;

ylabel('True Positive Rate');
xlabel('False Positive Rate')
title('MCI negative / MCI positive');

%Further post-processing the figure
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...
    'LineWidth'   , .5        );

axis square;

exportgraphics(f,config.ResultFolder+"/ROCL3_L3Hat_MCIposvsneg"+".png",'Resolution',300);
exportgraphics(f,config.ResultFolder+"/ROCL3_L3Hat_MCIposvsneg"+".pdf",'Resolution',300, 'ContentType','vector');

clear parametersName filesName i colors f ll legendText 

%% get the return distance and the actual distance for each trial
% each value is a mean value for per participant
% l3 is the real distance to the first cone (origin)
% l3hat is the actual distance of the last segme
function [l3_l3hat] = getReturnAndActualDistance(AllX, AllDX)
    
    nParticipants = length(AllX);

    %Averaged Across Trials
    l3      = nan(1,nParticipants);
    l3hat    = nan(1,nParticipants);

    for ip = 1:nParticipants

        X = AllX{ip};
        DX = AllDX{ip};
        
        currl3    = nan(1,length(X));
        currl3hat = nan(1,length(X));
        
        for tr = 1:length(X)
            currl3(tr)    = DX{tr}(3);
            currl3hat(tr) = norm(X{tr}(3,:));
        end

        l3(ip)    = mean(currl3);
        l3hat(ip) = mean(currl3hat);

    end

    l3_l3hat = l3' - l3hat';

end

%%
function AUC = plotROCCurve(param1, param2, paramName, filename, param1Label, param2Label, paramColor, config)
    % Preparing the logistic regression
    allData = [param1; param2];
    allDatalogicalResponse = (1:height(param1) + height(param2))' > height(param1);
    label1     = cell(height(param1),1);
    label1(:)  = {param1Label};
    label2     = cell(height(param2),1);
    label2(:)  = {param2Label};
    allLabels  = [label1;label2];

    clear hcLabel mciLabel
    
    % Fitting the logistic regression
    mdl = fitglm(allData,allDatalogicalResponse,'Distribution', 'binomial','Link','logit');
    
    allDataScores = mdl.Fitted.Probability;
    [X,Y,~,AUC.Value] = perfcurve(allLabels, allDataScores, param2Label, NBoot=10);

    plot(X(:,1),Y(:,1),...
        "Color",paramColor,...
        'LineWidth',1.0...
        );
end

%% A function for getting Results from All Conditions
function [AllParams, AllX, AllDX, AllTheta, AllIC] = getResultsAllConditions(TransformedData, config)
    %get the estimated parameters, X, DX, Theta, IC for all trial
    %conditions for each group of data.    
    config.TrialFilter = 0;
    tic
    disp('%%%%%%%%%%%%%%% PERFORMING FITTING %%%%%%%%%%%%%%%');
    Results = PerformGroupFit(TransformedData, config);
    toc

    AllParams = Results.estimatedParams;
    AllX      = Results.X;
    AllDX     = Results.DX;
    AllTheta  = Results.THETADX;
    AllIC     = Results.IC;
    
end