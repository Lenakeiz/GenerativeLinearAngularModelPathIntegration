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
%% Model fitting
%Model related parameters
config.ModelName        = "ConstSpeedModelDistAngle";
config.ParamName        = ["beta", "sigma", "g3", "nu"];
config.subtype          = "Regress2Mean"; %DistAng_SimpleGain, DistAng_RGmean

% config.ModelName        = "ConstSpeedModel5Params";
% config.ParamName        = ["beta", "g3", "b", "sigma", "nu"];
% config.subtype          = "DistAng_RGb"; %DistAng_SimpleGain, DistAng_RGmean

% config.ModelName        = "ConstSpeedModel";
% config.ParamName        = ["beta", "bG3", "g2", "g3", 'b', "sigma", "nu"];
% config.subtype          = "DistAng_RGmean";%choose from 1,egoNoise / 2, onlyDist / 3, onlyAng_RGb, 
%                                                       %4, onlyAng_RGmean / 5, DistAng_RGb / 6, DistAng_RGmean
config.includeStand     = false;
config.useOoBTrial      = false;
config.useweber         = false; %only true when use weber law in simple generative models
config.NumTotalParams   = length(config.ParamName);
config.NumFreeParams    = 4;

resultfolder = savefolder+"PaperFigs/OurDataResults/XFig5_"+config.ModelName+"_"+config.subtype+"_exclude2patcipants";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% Model fitting for MCIPos
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCIPos   = CalculateTrackingPath(MCIPos, config);
MCIPos = TransformPaths(MCIPos);%transform data
[AllMCIPosParams, AllMCIPosX, AllMCIPosDX, AllMCIPosTheta, AllMCIPosIC] = getResultsAllConditions(MCIPos, config);

%% Model fitting for MCINeg
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCINeg   = CalculateTrackingPath(MCINeg, config);
MCINeg = TransformPaths(MCINeg);%transform data
[AllMCINegParams,  ~, ~, ~, ~] = getResultsAllConditions(MCINeg, config);

%% Setting colors for using in plots
ColorPattern; 

%% TwowayAnova
[anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnova_LIModel_MCIPosvsMCINeg(AllMCIPosParams, AllMCINegParams, config);

%% BarScatter Plot between Young and HealthyOld for all Fitted Params
BoxPlotOfFittedParam(AllMCIPosParams, AllMCINegParams, anova_tab, config);
BoxPlotOfFittedParamMergeCondition(AllMCIPosParams, AllMCINegParams, multicomp_tab1, config)

%%
function BoxPlotOfFittedParam(AllMCIPosParams, AllMCINegParams, anova_tab, config)
    
    numConds = 3; %3 is the condition number
    ParamName = config.ParamName;
    for ParamIndx=1:length(ParamName)

        MCIPosParamAllConds = []; %dimension are different, so separate from MCIParamAllConds
        MCINegParamAllConds = [];
        for TRIAL_FILTER=1:numConds
            %% extract data
            MCIPosParam = AllMCIPosParams{TRIAL_FILTER}(:,ParamIndx);
            MCIPosParamAllConds = [MCIPosParamAllConds,MCIPosParam];

            MCINegParam = AllMCINegParams{TRIAL_FILTER}(:,ParamIndx);
            MCINegParamAllConds = [MCINegParamAllConds,MCINegParam];
        end

        %remove Nan rows (nan coz of 1, removing participants with short walking length; 2, not enough trials for parameter estimation)
        MCIPosParamAllConds = removeNanRows(MCIPosParamAllConds);
        MCINegParamAllConds = removeNanRows(MCINegParamAllConds);
    
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
        colorForMCIPos = config.color_scheme_npg(6,:);
        colorForMCINeg = config.color_scheme_npg(3,:);

        %set params
        whisker_value               =   1.5;
        box_lineWidth               =   0.3;
        box_widths_value            =   0.2;
        box_color_transparency      =   0.5;        %faceAlpha
        center_x                    =   [1,2,3];    %center of box (three conditions)
        shift_value                 =   0.2;        %box shift from center
        median_lineWidth            =   2;
        median_color                =   'k';
        scatter_jitter_value        =   0.1;
        scatter_markerSize          =   10;
        scatter_marker_edgeColor    =   'k';
        scatter_marker_edgeWidth    =   0.5;
        scatter_color_transparency  =   0.7;        %faceAlpha        

        %% boxplot for each column in MCIPos
        bp1 = boxplot(MCIPosParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x-shift_value);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %boxplot for each column in MCIMerge
        bp2 = boxplot(MCINegParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x+shift_value);
        set(bp2,'linewidth',box_lineWidth);

        %% Coloring each box
        %findobj first getting the three boxes for MCINeg (from bp2) from right to left
        %then getting the three boxes for MCIPos (frm bp1) from right to left
        h = findobj(gca,'Tag','Box'); 
        for i = 1:length(h)
            if i<4  %get the MCINeg box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForMCINeg,'FaceAlpha',box_color_transparency);
            else %get the MCIPos box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForMCIPos,'FaceAlpha',box_color_transparency);
            end
        end

        %% Adjusting median
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% add scatter plot and the mean of MCIPos
        num_points = size(MCIPosParamAllConds,1);
        for i=1:size(MCIPosParamAllConds,2)
            hold on
            x = i*ones(num_points,1)-shift_value+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, MCIPosParamAllConds(:,i), scatter_markerSize, ...
                    'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForMCIPos, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
            %add errorbar
            mean_MCIPos = mean(MCIPosParamAllConds(:,i));
            sem_MCIPos = std(MCIPosParamAllConds(:,i))./sqrt(length(MCIPosParamAllConds(:,i)));
            errorbar(i-shift_value,mean_MCIPos,sem_MCIPos,'k','LineStyle','None', 'LineWidth', 2,  'CapSize', 14);    
            hold on
            %add mean point
            scatter(i-shift_value, mean_MCIPos, 4*scatter_markerSize, 'd',...
                    'filled','MarkerEdgeColor','k', ...
                    'MarkerFaceColor','w', ...
                    'LineWidth',scatter_marker_edgeWidth);
        end

        %% add scatter plot and the mean of MCINeg
        num_points = size(MCINegParamAllConds,1);
        for i=1:size(MCINegParamAllConds,2)
            hold on
            x = i*ones(num_points,1)+shift_value+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, MCINegParamAllConds(:,i), scatter_markerSize, ...
                    'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForMCINeg, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
            %add errorbar
            mean_MCI = mean(MCINegParamAllConds(:,i));
            sem_MCI = std(MCINegParamAllConds(:,i))./sqrt(length(MCINegParamAllConds(:,i)));
            errorbar(i+shift_value,mean_MCI,sem_MCI,'k','LineStyle','None', 'LineWidth', 2,  'CapSize', 14); 
            hold on
            %add mean point
            scatter(i+shift_value, mean_MCI, 4*scatter_markerSize, 'd',...
                    'filled','MarkerEdgeColor','k', ...
                    'MarkerFaceColor','w', ...
                    'LineWidth',scatter_marker_edgeWidth);
        end

        %% Further post-processing the figure

        %calculate te ylim
        alldata = [MCIPosParamAllConds;MCINegParamAllConds];
        maxdata = max(alldata,[],'all');
        mindata = min(alldata, [], 'all');
        lowupYlim = [mindata-.1*(maxdata-mindata)-eps, maxdata+.1*(maxdata-mindata)+eps]; 
        %eps to make sure when mindata=maxdata, there won't be an error

        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.1 .1 .1], ...
            'YColor'      , [.1 .1 .1], ...
            'XTick'       , (1:3),... 
            'XLim'        , [0.5, 3.5],...
            'YLim'        , lowupYlim,...   
            'XTickLabel'  , {'No Change','No Distal Cue', 'No Optical Flow'},...
            'LineWidth'   , .5        );
        ylabel(ParamName(ParamIndx));
        allpatches = findall(gca,'type','Patch');
        legend(allpatches(1:3:end), {'MCIPos' 'MCINeg'}, 'Location','northeast', 'NumColumns',2);

        %extract pvalue for group, conditino and interaction to show on the figure 
        anova_result = anova_tab{ParamIndx};
        group_pvalue = anova_result{2,7};
        condition_pvalue = anova_result{3,7};
        interaction_pvalue = anova_result{4,7};
%         str = {['Group Pvalue = ',sprintf('%.2g',group_pvalue)],...
%                ['Condition Pvalue = ',sprintf('%.2g',condition_pvalue)],...
%                ['Interaction Pvalue = ',sprintf('%.2g',interaction_pvalue)]};
%         annotation('textbox',[0.2 0.6 0.3 0.3],'String',str,'FitBoxToText','on');
        title(strcat(['Group P = ',sprintf('%.2g',group_pvalue)],...
              ['    Condition P = ',sprintf('%.2g',condition_pvalue)],...
              ['    Interaction P = ',sprintf('%.2g',interaction_pvalue)]))

        %% save figure
        exportgraphics(f,config.ResultFolder+"/Box_"+ParamName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/Box_"+ParamName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end

%%
function BoxPlotOfFittedParamMergeCondition(AllMCIPosParams, AllMCINegParams, multicomp_tab1, config)
    
    numConds = 3; %3 is the condition number
    ParamName = config.ParamName;
    for ParamIndx=1:length(ParamName)

        MCIPosParamAllConds = []; 
        MCINegParamAllConds = [];

        for TRIAL_FILTER=1:numConds
            %% extract data
            MCIPosParam = AllMCIPosParams{TRIAL_FILTER}(:,ParamIndx);
            MCIPosParamAllConds = [MCIPosParamAllConds,MCIPosParam];

            MCINegParam = AllMCINegParams{TRIAL_FILTER}(:,ParamIndx);
            MCINegParamAllConds = [MCINegParamAllConds,MCINegParam];
        end

        %remove Nan rows (nan coz of 1, removing participants with short walking length; 2, not enough trials for parameter estimation)
        MCIPosParamAllConds = removeNanRows(MCIPosParamAllConds);
        MCINegParamAllConds = removeNanRows(MCINegParamAllConds);

        MCIPosParamMean = mean(MCIPosParamAllConds, 2);
        MCINegParamMean = mean(MCINegParamAllConds, 2);
    
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
        whisker_value               =   1.5;
        box_lineWidth               =   0.3;
        box_widths_value            =   0.4;
        box_color_transparency      =   0.5; %faceAlpha
        median_lineWidth            =   2;
        median_color                =   'k';
        scatter_jitter_value        =   0.2;
        scatter_markerSize          =   10;
        scatter_marker_edgeColor    =   'k';
        scatter_marker_edgeWidth    =   0.5;
        scatter_color_transparency  =   0.7; %faceAlpha        

        %% boxplot for each column in MCIPOs
        bp1 = boxplot(MCIPosParamMean, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 1);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %boxplot for each column in MCINeg
        bp2 = boxplot(MCINegParamMean, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 2);
        set(bp2,'linewidth',box_lineWidth);

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
        num_points = length(MCIPosParamMean);
        hold on
        x = ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, MCIPosParamMean, scatter_markerSize, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForMCIPos, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        %add errorbar
        mean_MCIPos = mean(MCIPosParamMean);
        sem_MCIPos= std(MCIPosParamMean)./sqrt(length(MCIPosParamMean));
        errorbar(1,mean_MCIPos,sem_MCIPos,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18);    
        hold on
        %add mean point
        scatter(1, mean_MCIPos, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);

        %% add scatter plot and the mean of MCI
        num_points = size(MCINegParamMean,1);
        hold on
        x = 2*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, MCINegParamMean, scatter_markerSize, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForMCINeg, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        %add errorbar
        mean_MCINeg = mean(MCINegParamMean);
        sem_MCINeg = std(MCINegParamMean)./sqrt(length(MCINegParamMean));
        errorbar(2,mean_MCINeg,sem_MCINeg,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18); 
        hold on
        %add mean point
        scatter(2, mean_MCINeg, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);        

        %% Further post-processing the figure
        %calculate the ylim
        alldata = [MCIPosParamMean;MCINegParamMean];
        maxdata = max(alldata,[],'all');
        mindata = min(alldata, [], 'all');
        lowupYlim = [mindata-.1*(maxdata-mindata)-eps, maxdata+.1*(maxdata-mindata)+eps]; 
        %eps to make sure when mindata=maxdata, there won't be an error

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
        ylabel(ParamName(ParamIndx));

        %extract pvalue for multicomparison of Group effect for showing on the figure
        multicomp_result = multicomp_tab1{ParamIndx};
        Pvalue = multicomp_result(1,6); % MCIPos vs. MCINeg see Two-way anova for details
%         str = {['P = ',sprintf('%.2g',Pvalue)]};
%         annotation('textbox',[0.2 0.6 0.3 0.3],'String',str,'FitBoxToText','on');
        title(['P value = ',sprintf('%.2g',Pvalue)])

        %% add sigstar 
        if Pvalue<0.05
            H=sigstar({[1,2]},[Pvalue]);
        end

        %% save figure
        exportgraphics(f,config.ResultFolder+"/ZMergeCondsBox_"+ParamName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/ZMergeCondsBox_"+ParamName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end

