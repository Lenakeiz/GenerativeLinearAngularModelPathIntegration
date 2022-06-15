%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default'); %for code reproducibility

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
savefolder = pwd + "/Output/";

%% setting the configuration
config.Speed.alpha                                  = 0.9;                 %Paramanter for running speed calculation
config.Speed.timeOffsetAfterFlagReach               = 1.5;                 %Time to track after flag reached in seconds 
config.Speed.smoothWindow                           = 10;                  % tracking rate should be 10Hz so 4 secs window is 40 datapoints
config.Speed.velocityCutoff                         = 0.2;                 % velocity cutoff to select only the walking part of the reconstructed velocity
config.Speed.timeOffsetForDetectedTemporalWindow    = 0.4;                 % time in seconds that will push earlier/ the detected rising edge
config.UseGlobalSearch                              = true;
config.TrackedInboundAngularDeltaT                  = 1;
config.TrialFilter                                  = 0;                   %merge three conditions

%% Model fitting

% config.ModelName        = "ConstSpeedModelwith_g2";
% config.ParamName        = ["beta", "g2", "g3", "sigma", "nu"];

config.ModelName        = "ConstSpeedModelwith_g2_k3";
config.ParamName        = ["beta", "g2", "g3", "k3", "sigma", "nu"];

config.includeStand     = false;
config.useweber         = false; %only true when use weber law in simple generative models
config.NumParams        = length(config.ParamName);

resultfolder = savefolder+"PaperFigs/ModelAfterDataCleaning/Fig3_"+config.ModelName;
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% Model fitting for YoungControl data
config.Speed.tresholdForBadParticipantL1Recontruction = 1.55;   % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
%transform data
YoungControls = TransformPaths(YoungControls);
YoungControls   = CalculateTrackingPath(YoungControls, config);
ManuallyScoringYoung;
YoungResults  = PerformGroupFit(YoungControls, config);
YoungParams = YoungResults.estimatedParams;

%% Model fitting for HealthyOld data
config.Speed.tresholdForBadParticipantL1Recontruction = 2.0; 
%transform data
HealthyControls = TransformPaths(HealthyControls);
HealthyControls   = CalculateTrackingPath(HealthyControls, config);
ManuallyScoringHealthyOld;
HealthyOldResults  = PerformGroupFit(HealthyControls, config);
HealthyOldParams = HealthyOldResults.estimatedParams;

%% Model fitting for MCIPos
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
%transform data
MCIPos = TransformPaths(MCIPos);
MCIPos   = CalculateTrackingPath(MCIPos, config);
ManuallyScoringMCIPos;
MCIPosResults  = PerformGroupFit(MCIPos, config);
MCIPosParams = MCIPosResults.estimatedParams;

%% Model fitting for MCINeg
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
%transform data
MCINeg = TransformPaths(MCINeg);
MCINeg   = CalculateTrackingPath(MCINeg, config);
ManuallyScoringMCINeg;
MCINegResults  = PerformGroupFit(MCINeg, config);
MCINegParams = MCINegResults.estimatedParams;

%% Model fitting for MCIUnk
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
%transform data
MCIUnk = TransformPaths(Unknown);
MCIUnk   = CalculateTrackingPath(MCIUnk, config);
ManuallyScoringMCIUnk;
MCIUnkResults  = PerformGroupFit(MCIUnk, config);
MCIUnkParams = MCIUnkResults.estimatedParams;

%% Setting colors for using in plots
ColorPattern; 

%% merge MCI together
MCIParams = [MCIPosParams;MCINegParams;MCIUnkParams];

%% OnewayAnova
[~,multicomp_tab] = OnewayAnova_Young_HealthyOld_MergeMCI(YoungParams, HealthyOldParams, MCIParams, config);

%% BarScatter Plot between Young and HealthyOld for all Fitted Params
plotBoxOfFittedParam(YoungParams, HealthyOldParams, MCIParams, multicomp_tab, config);

%%
function plotBoxOfFittedParam(YoungParams, HealthyOldParams, MCIParams, multicomp_tab, config)
    
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
        colorForYoung = config.color_scheme_npg(3,:);        
        colorForHOld = config.color_scheme_npg(5,:);
        colorForMCI = config.color_scheme_npg(2,:);

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

        %% boxplot for each column in MCIMerge
        bp1 = boxplot(MCIParams(:, ParamIndx), ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 3);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %% boxplot for each column in HealthOld
        bp2 = boxplot(HealthyOldParams(:,ParamIndx), ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 2);
        set(bp2,'linewidth',box_lineWidth);

        hold on
        %% boxplot for each column in Young
        bp3 = boxplot(YoungParams(:,ParamIndx), ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 1);
        set(bp3,'linewidth',box_lineWidth);
        

        %% Coloring each box
        %findobj first getting the box for Young (from bp3) 
        % then for HealthyOld (from bp2) 
        %then getting the box for MCI (frm bp1)
        h = findobj(gca,'Tag','Box'); 
        %get the Young box
        patch(get(h(1),'XData'),get(h(1),'YData'),colorForYoung,'FaceAlpha',box_color_transparency);        
        %get the HelthyOld box
        patch(get(h(2),'XData'),get(h(2),'YData'),colorForHOld,'FaceAlpha',box_color_transparency);
        %get the MCI box
        patch(get(h(3),'XData'),get(h(3),'YData'),colorForMCI,'FaceAlpha',box_color_transparency);

        %% Adjusting median
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% add scatter plot and the mean of Young
        num_points = size(YoungParams,1);
        x = 1*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        hold on
        scatter(x, YoungParams(:,ParamIndx), scatter_markerSize, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForYoung, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        %add errorbar
        mean_Young = mean(YoungParams(:,ParamIndx));
        sem_Young = std(YoungParams(:,ParamIndx))./sqrt(length(YoungParams(:,ParamIndx)));
        errorbar(1,mean_Young,sem_Young,'k','LineStyle','None', 'LineWidth', 2); 
        hold on
        %add mean point
        scatter(1, mean_Young, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);    

        %% add scatter plot and the mean of HealthyOld
        num_points = size(HealthyOldParams,1);
        hold on
        x = 2*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, HealthyOldParams(:,ParamIndx), scatter_markerSize, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForHOld, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        %add errorbar
        mean_Hold = mean(HealthyOldParams(:,ParamIndx));
        sem_Hold = std(HealthyOldParams(:,ParamIndx))./sqrt(length(HealthyOldParams(:,ParamIndx)));
        errorbar(2,mean_Hold,sem_Hold,'k','LineStyle','None', 'LineWidth', 2);    
        hold on
        %add mean point
        scatter(2, mean_Hold, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);    

        %% add scatter plot and the mean of MCI
        num_points = size(MCIParams,1);
        hold on
        x = 3*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, MCIParams(:,ParamIndx), scatter_markerSize, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForMCI, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        %add errorbar
        mean_MCI = mean(MCIParams(:,ParamIndx));
        sem_MCI = std(MCIParams(:,ParamIndx))./sqrt(length(MCIParams(:,ParamIndx)));
        errorbar(3,mean_MCI,sem_MCI,'k','LineStyle','None', 'LineWidth', 2); 
        hold on
        %add mean point
        scatter(3, mean_MCI, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);     

        %% Further post-processing the figure

        %calculate the ylim
        alldata = [YoungParams;HealthyOldParams;MCIParams];
        maxdata = max(alldata,[],'all');
        mindata = min(alldata, [], 'all');
        lowupYlim = [mindata-.1*(maxdata-mindata)-eps, maxdata+.1*(maxdata-mindata)+eps]; 

        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.1 .1 .1], ...
            'YColor'      , [.1 .1 .1], ...
            'XTick'       , (1:3),... 
            'XLim'        , [0.5, 3.5],...
            'YLim'        , lowupYlim,...
            'XTickLabel'  , {'Young','HealthyOld','MCI'},...
            'LineWidth'   , .5        );
            %'Ytick'       , [0,0.5,1.0,1.5],...
            %'YLim'        , [0, 1.5],...   
        ylabel(ParamName(ParamIndx));
        legend(gca, {'Young', 'HealthyOld', 'MCI'}, 'Location','northeast', 'NumColumns',1);
        %xlabel('G_1','Interpreter','tex'); ylabel('G_2','Interpreter','tex');

        %extract pvalue for multicomparison of Group effect for showing on the figure
        multicomp_result = multicomp_tab{ParamIndx};
        % 1       2             3
        % MCI     HealthyOld    Young
        PvalueYoungvsHealthyOld = multicomp_result(3,6); % Young vs. HealthyOld see Two-way anova for details
        PvalueHealthyOldvsMCI = multicomp_result(1,6); % HealthyOld v.s. MCI vs.  see Two-way anova for details
        PvalueYoungvsMCI = multicomp_result(2,6); % Young vs. MCI see Two-way anova for details 
        title(strcat(['P12 = ',sprintf('%.2g',PvalueYoungvsHealthyOld)],...
              ['    P23 = ',sprintf('%.2g',PvalueHealthyOldvsMCI)],...
              ['    P13 = ',sprintf('%.2g',PvalueYoungvsMCI)]))

        %% add sigstar 
        AllP = [PvalueYoungvsHealthyOld,PvalueHealthyOldvsMCI,PvalueYoungvsMCI];
        Xval = [[1,2];[2,3];[1,3]];
        %select those P value smaller than 0.05 (only add line when p<0.05)
        % * represents p<=0.05
        % ** represents p<=1E-2
        % *** represents p<=1E-3
        PsigInd = AllP<0.05;
        if sum(AllP<0.05)>0
            Xval_select = Xval(PsigInd,:);
            AllP_select = AllP(PsigInd);
            XXX = {};
            for i=1:sum(AllP<0.05)
                XXX{i} = Xval_select(i,:);
            end
            H=sigstar(XXX,AllP_select);
        end

        %% save figure
        exportgraphics(f,config.ResultFolder+"/ZMergeCondsBox_"+ParamName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/ZMergeCondsBox_"+ParamName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end
