%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default'); %for code reproducibility

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
savefolder = pwd + "/Output/";

%% setting the configuration
config.Speed.alpha                                  = 0.9;                 % Parameter for running speed calculation
config.Speed.timeOffsetAfterFlagReach               = 1.5;                 % Time to track after flag reached in seconds 
config.Speed.smoothWindow                           = 10;                  % tracking rate should be 10Hz so 4 secs window is 40 datapoints
config.Speed.velocityCutoff                         = 0.2;                 % velocity cutoff to select only the walking part of the reconstructed velocity
config.Speed.timeOffsetForDetectedTemporalWindow    = 0.4;                 % time in seconds that will push earlier/ the detected rising edge
config.UseGlobalSearch                              = true;
config.TrackedInboundAngularDeltaT                  = 1;
config.TrialFilter                                  = 0;                   %merge three conditions

%% Model fitting

config.ModelName        = "ConstSpeedModelwith_g2";
config.ParamName        = ["beta", "g2", "g3", "sigma", "nu"];

% config.ModelName        = "ConstSpeedModelwith_g2_k3";
% config.ParamName        = ["beta", "g2", "g3", "k3", "sigma", "nu"];

config.includeStand     = false;
config.useweber         = false; %only true when use weber law in simple generative models
config.NumParams   = length(config.ParamName);

resultfolder = savefolder+"PaperFigs/ModelAfterDataCleaning/Fig4_"+config.ModelName;
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% Model fitting for MCIPos
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCIPos = TransformPaths(MCIPos);%transform data
MCIPos   = CalculateTrackingPath(MCIPos, config);
ManuallyScoringMCIPos;
MCIPosResults  = PerformGroupFit(MCIPos, config);
MCIPosParams = MCIPosResults.estimatedParams;

%% Model fitting for MCINeg
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0;
MCINeg = TransformPaths(MCINeg);%transform data
MCINeg   = CalculateTrackingPath(MCINeg, config);
ManuallyScoringMCINeg;
MCINegResults  = PerformGroupFit(MCINeg, config);
MCINegParams = MCINegResults.estimatedParams;

%% Setting colors for using in plots
ColorPattern; 

%% OnewayAnova
[~,multicomp_tab] = OnewayAnova_MCIPosvsMCINeg(MCIPosParams, MCINegParams, config);

%% BarScatter Plot 
plotBoxOfFittedParam(MCIPosParams, MCINegParams, multicomp_tab, config);

%%
function plotBoxOfFittedParam(MCIPosParams, MCINegParams, multicomp_tab, config)
    
    ParamName = config.ParamName;

    for ParamIndx=1:length(ParamName)

        %% set figure info
        f = figure('visible','off','Position', [100 100 500 500]);
        %%% Font type and size setting %%%
        % Using Arial as default because all journals normally require the font to
        % be either Arial or Helvetica
        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)     

        %%% Color definition %%%
        colorForMCIPos = config.color_scheme_npg(6,:);
        colorForMCINeg = config.color_scheme_npg(3,:);

        %set params
        whisker_value = 1.5;
        box_lineWidth = 0.3;
        box_widths_value = 0.2;
        box_color_transparency = 0.5; %faceAlpha
        median_lineWidth = 2;
        median_color = 'k';
        scatter_jitter_value = 0.1;
        scatter_markerSize=10;
        scatter_marker_edgeColor = 'k';
        scatter_marker_edgeWidth = 0.5;
        scatter_color_transparency = 0.7; %faceAlpha        

        %% boxplot for each column in MCINeg
        bp1 = boxplot(MCINegParams(:, ParamIndx), ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 2);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %% boxplot for each column in MCIPos
        bp2 = boxplot(MCIPosParams(:,ParamIndx), ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 1);
        set(bp2,'linewidth',box_lineWidth);

        hold on       

        %% Coloring each box
        %findobj first getting the box for MCINeg (from bp2) 
        %then getting the box for MCIPos(frm bp1)
        h = findobj(gca,'Tag','Box'); 
        %get the MCI box
        patch(get(h(1),'XData'),get(h(1),'YData'),colorForMCINeg,'FaceAlpha',box_color_transparency);
        %get the HelthyOld box
        patch(get(h(2),'XData'),get(h(2),'YData'),colorForMCIPos,'FaceAlpha',box_color_transparency);

        %% Adjusting median
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% add scatter plot and the mean of MCIPos
        num_points = size(MCIPosParams,1);
        hold on
        x = 1*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, MCIPosParams(:,ParamIndx), scatter_markerSize, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForMCIPos, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        %add errorbar
        mean_MCIPos = mean(MCIPosParams(:,ParamIndx));
        sem_MCIPos = std(MCIPosParams(:,ParamIndx))./sqrt(length(MCIPosParams(:,ParamIndx)));
        errorbar(1,mean_MCIPos,sem_MCIPos,'k','LineStyle','None', 'LineWidth', 2);    
        hold on
        %add mean point
        scatter(1, mean_MCIPos, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);

        %% add scatter plot and the mean of MCI
        num_points = size(MCINegParams,1);
        hold on
        x = 2*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, MCINegParams(:, ParamIndx), scatter_markerSize, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForMCINeg, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        %add errorbar
        mean_MCINeg = mean(MCINegParams(:, ParamIndx));
        sem_MCINeg = std(MCINegParams(:, ParamIndx))./sqrt(length(MCINegParams(:, ParamIndx)));
        errorbar(2,mean_MCINeg,sem_MCINeg,'k','LineStyle','None', 'LineWidth', 2); 
        hold on
        %add mean point
        scatter(2, mean_MCINeg, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);     

        %% Further post-processing the figure
        alldata = [MCIPosParams(:,ParamIndx); MCINegParams(:,ParamIndx)];
        maxdata = max(alldata);
        mindata = min(alldata);

        lowupYlim = [mindata-.2*(maxdata-mindata)-eps, maxdata+.4*(maxdata-mindata)+eps]; 

        %% Further post-processing the figure
        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.1 .1 .1], ...
            'YColor'      , [.1 .1 .1], ...
            'XTick'       , (1:2),... 
            'XLim'        , [0.5, 2.5],...
            'YLim'        , lowupYlim,...
            'XTickLabel'  , {'MCIPos','MCINeg'},...
            'LineWidth'   , .5        );
            %'Ytick'       , [0,0.5,1.0,1.5],... 
        ylabel(ParamName(ParamIndx));

        %% extract pvalue for multicomparison of Group effect for showing on the figure
        multicomp_result = multicomp_tab{ParamIndx};
        Pvalue = multicomp_result(1,6); %POsthoc test Pvalue between MCIPos and MCINeg 
        %str = {['Pvalue = ',sprintf('%.2e',Pvalue)]};
        %annotation('textbox',[0.2 0.6 0.3 0.3],'String',str,'FitBoxToText','on');
        title(['Pvalue = ',sprintf('%.2e',Pvalue)])
        hold on

        if Pvalue<0.05
            % add sigstar
            H=sigstar({[1,2]},[Pvalue]);
        end

        %% save figure
        exportgraphics(f,config.ResultFolder+"/ZMergeCondsBox_"+ParamName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/ZMergeCondsBox_"+ParamName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end

