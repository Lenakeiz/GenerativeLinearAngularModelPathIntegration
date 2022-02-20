%% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
YoungControls = TransformPaths(YoungControls);
savefolder = pwd + "/Output/";

%% setting the configuration
config.UseGlobalSearch = true;
resultfolder = savefolder+"PaperFigs/Fig3";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% Model fitting using G1G2 model with only Distance error
%Model parameter G1 G2, sigma, nu. #params=4
config.ModelName = "DistErrG1G2";
config.NumParams = 4;
[AllYoungParams, AllYoungX, AllYoungDX, AllYoungTheta, AllYoungIC] = getResultsAllConditions(YoungControls, config);

%% Setting colors for using in plots
ColorPattern; 

%% BarScatter Plot between G1 and G2
plotBarScatterBetweenG1G2(AllYoungParams, config)

%% Correlation plot between G1 and G2
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

function plotBarScatterBetweenG1G2(AllParams, config)
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
        color_scheme_npg = config.color_scheme_npg(2,:);
        pairline_color = 'k';

        %% extract data
        ParamValues = AllParams{TRIAL_FILTER};
        subjectSize = size(ParamValues,1);
        all_G1 = zeros(1,subjectSize); all_G2 = zeros(1,subjectSize);
        for subj=1:subjectSize %for each subject
            paramX = ParamValues(subj,:);
            %extract G1 and G2
            G1 = paramX(1); all_G1(subj)=G1;
            G2 = paramX(2); all_G2(subj)=G2;

            jitter_value = 0.2*(rand(1)-0.5);
            %link one pair of points horizontally
            pt = plot([1+jitter_value, 2+jitter_value], [G1, G2], '-', 'Color', pairline_color ,'linewidth',2);            
            %transparancy
            pt.Color = [pt.Color 0.1];
            hold on
            %scatter
            st = scatter([1+jitter_value, 2+jitter_value], [G1, G2], 50, ...
                "filled", 'MarkerEdgeColor', 'k', 'MarkerFaceColor', pairline_color, ...
                'MarkerEdgeAlpha', 0.8, 'MarkerFaceAlpha', 0.3);
            hold on
        end
        meanG1 = mean(all_G1); meanG2 = mean(all_G2);
        b = bar([1,2],[meanG1,meanG2],'BarWidth', 0.4, ...
            'FaceColor', color_scheme_npg, 'FaceAlpha',0.5, 'LineWidth',1);

        hold on
        % Plot the errorbars
        stdG1 = std(all_G1);
        stdG2 = std(all_G2);
        eb = errorbar([1,2],[meanG1,meanG2],[stdG1,stdG2],'k', 'LineStyle', 'None', 'LineWidth', 2);

        % Sending barplot and errorbar plot back of the figure so that scatters can be drawn on top of it
        set(gca,'children',flipud(get(gca,'children'))) 

        %link the mean value
        pt = plot([1,2],[meanG1, meanG2], 'd-', 'Color', pairline_color ,'linewidth',6);
        pt.Color = [pt.Color 0.9];

        hold on
        %add y=1 for reference
        yline(1, 'LineStyle','-.', 'LineWidth', 1, 'Color', config.color_scheme_npg(1,:));

        %Further post-processing the figure
        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.1 .1 .1], ...
            'YColor'      , [.1 .1 .1], ...
            'XLim'        , [0.5, 2.5],...
            'YLim'        , [0, 1.5],...   
            'XTick'       , [1:2],... 
            'XTickLabel'  , {'G_1','G_2'},...
            'Ytick'       , [0,0.5,1.0,1.5],...
            'LineWidth'   , .5        );
        ylabel('Estimated value');
        %xlabel('G_1','Interpreter','tex'); ylabel('G_2','Interpreter','tex');
        
        %% save figure
        exportgraphics(f,config.ResultFolder+"/BarG1G2"+conditionName{TRIAL_FILTER}+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/BarG1G2"+conditionName{TRIAL_FILTER}+".pdf",'Resolution',300, 'ContentType','vector');
    end
end

function plotRegressionBetweenG1G2Square(AllParams, config)
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
        scattercolor = config.color_scheme_npg(2,:);

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
        set(h(1), 'Color', 'w', 'Marker', '.', 'MarkerSize', 0.1); %set data point making it invisible reset in scatter plot below
        set(h(2), 'Color', 'k', 'LineWidth', 3); %set the regression line
        set(h(3), 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.'); %set the upper confidence line
        set(h(4), 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.'); %set the lower confidence line
        
        hold on
        %patch between confidence lines
        patch([h(3).XData,fliplr(h(4).XData)],[h(3).YData, fliplr(h(4).YData)], ...
               'k', 'edgecolor', 'none', 'facealpha', 0.1);
        %refline(1,0); 

        hold on 
        scatter(all_G1, all_G2.^2, 50, "filled", 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', scattercolor, 'MarkerEdgeAlpha', 1, 'MarkerFaceAlpha', 0.5)
        
        [R,P]=corrcoef(all_G1,all_G2.^2);
        str = {['Person Corrcoef = ',num2str(round(R(1,2),4))],...
               ['Pvalue = ',num2str(round(P(1,2),4))]};
        annotation('textbox',[0.2 0.6 0.3 0.3],'String',str,'FitBoxToText','on');

        %Further post-processing the figure
        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.1 .1 .1], ...
            'YColor'      , [.1 .1 .1], ...
            'XLim'        , [0, 1.5],...
            'YLim'        , [0, 1.5],...   
            'XTick'       , [0,0.5,1.0,1.5],... 
            'Ytick'       , [0,0.5,1.0,1.5],...            
            'LineWidth'   , .5        );

        xlabel('G_1','Interpreter','tex'); ylabel('G_2^2','Interpreter','tex');
        legend('off');
        
        title(conditionName{TRIAL_FILTER});
        %% save figure
        exportgraphics(f,config.ResultFolder+"/G1_G2Square_"+conditionName{TRIAL_FILTER}+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/G1_G2Square_"+conditionName{TRIAL_FILTER}+".pdf",'Resolution',300, 'ContentType','vector');
    end
end



