%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default'); %for code reproducibility

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');
savefolder = pwd + "/Output/";

%% setting the configuration
config.Speed.alpha = 0.9;                                       %Paramanter for running speed calculation
config.Speed.timeOffsetAfterFlagReach = 2;                      %Time to track after flag reached in seconds 
config.Speed.smoothWindow = 10;                                 % tracking rate should be 10Hz so 4 secs window is 40 datapoints
config.Speed.velocityCutoff = 0.2;                              % velocity cutoff to select only the walking part of the reconstructed velocity
config.Speed.timeOffsetForDetectedTemporalWindow = 0.5;         % time in seconds that will push earlier/ the detected rising edge

config.TrialFilter = 0; %merge all conditions
config.UseGlobalSearch = true;

resultfolder = savefolder+"PaperFigs/Fig4";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end


%% Model fitting using Gamma Base Model
%Model parameter gamma, g3, b, sigma, nu. #params=5
config.ModelName = "LIFull";
config.NumParams = 5;

%% Model fitting for YoungControl data
%% calculating tracking path and transoform data
config.Speed.tresholdForBadParticipantL1Recontruction = 1.55;   % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
YoungControls   = CalculateTrackingPath(YoungControls, config);
%transform data
YoungControls = TransformPaths(YoungControls);
YoungmeanS = getMeanWalkingSpeed(YoungControls);

%% Model fitting for HealthyOld data
config.Speed.tresholdForBadParticipantL1Recontruction = 2.0; 
HealthyControls   = CalculateTrackingPath(HealthyControls, config);
%transform data
HealthyOld = TransformPaths(HealthyControls);
HealthyOldmeanS = getMeanWalkingSpeed(HealthyOld);

%% Model fitting for MCIPos
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCIPos   = CalculateTrackingPath(MCIPos, config);
%transform data
MCIPos = TransformPaths(MCIPos);
MCIPosmeanS = getMeanWalkingSpeed(MCIPos);

%% Model fitting for MCINeg
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCINeg   = CalculateTrackingPath(MCINeg, config);
%transform data
MCINeg = TransformPaths(MCINeg);
MCINegmeanS = getMeanWalkingSpeed(MCINeg);

%% Model fitting for MCIUnk
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
Unknown   = CalculateTrackingPath(Unknown, config);
%transform data
MCIUnk = TransformPaths(Unknown);
MCIUnkmeanS = getMeanWalkingSpeed(MCIUnk);

%% Setting colors for using in plots
ColorPattern; 

%% TwowayAnova
[anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnovaOn_Young_HealthyOld_MergeMCI(AllYoungParams, AllHealthyOldParams, AllMCIParams, config);

%% BarScatter Plot between Young and HealthyOld for all Fitted Params
plotBoxOfFittedParam(AllYoungParams, AllHealthyOldParams, AllMCIParams, anova_tab, config);
plotBoxOfFittedParamMergeCondition(AllYoungParams, AllHealthyOldParams, AllMCIParams, multicomp_tab1, config)

%% getting Mean Walking Speed of all participants in a group
function meanS = getMeanWalkingSpeed(Dat)
    numSubjs = length(Dat.Reconstructed)
    meanS = [];
    for subj=1:numSubjs
        %length on leg1
        length1 = table2array(Dat{subj}(:,'L1Real'));
        %duration on leg1
        duration1 = table2array(Dat{subj}(:,'T_L1'));
        meanspeed1 = length1./duration1;

        %length on leg2
        length2 = table2array(Dat{subj}(:,'L2Real'));
        %duration on leg2
        duration2 = table2array(Dat{subj}(:,'T_L2'));
        meanspeed2 = length2./duration2;
        
        %mean speed across trials of one participant
        meanspeed = mean([meanspeed1;meanspeed2]);
        
        meanS = [meanS,meanspeed];
    end
end

%%
function plotBoxOfFittedParam(AllYoungParams, AllHealthyOldParams, AllMCIParams, anova_tab, config)
    
    numConds = 3; %3 is the condition number
    ParamName = ["beta", "bG_3", "g_2", "Angular gain g_3", "Anglar bias b", "sigma", "nu"];
    StoreName = ["beta", "bG_3", "g_2", "g_3", "b", "sigma", "nu"];
    for ParamIndx=1:length(ParamName)

        YoungParamAllConds = [];
        HealthyOldParamAllConds = []; %dimension are different, so separate from YoungParams
        MCIParamAllConds = [];
        
        for TRIAL_FILTER=1:numConds
            %% extract data
            YoungParam = AllYoungParams{TRIAL_FILTER}(:,ParamIndx);
            YoungParamAllConds = [YoungParamAllConds,YoungParam];

            HealthyOldParam = AllHealthyOldParams{TRIAL_FILTER}(:,ParamIndx);
            HealthyOldParamAllConds = [HealthyOldParamAllConds,HealthyOldParam];

            MCIParam = AllMCIParams{TRIAL_FILTER}(:,ParamIndx);
            MCIParamAllConds = [MCIParamAllConds,MCIParam];            
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
        colorForYoung = config.color_scheme_npg(3,:);
        colorForHOld = config.color_scheme_npg(5,:);
        colorForMCI = config.color_scheme_npg(2,:);
        
        %set params
        whisker_value = 1.5;
        box_lineWidth = 0.3;
        box_widths_value = 0.2;
        box_color_transparency = 0.5; %faceAlpha
        %center of box (three conditions)
        center_x = [1,2,3];
        shift_value = 0.25; %box shift from center
        median_lineWidth = 2;
        median_color = 'k';
        scatter_jitter_value = 0.1;
        scatter_markerSize=10;
        scatter_marker_edgeColor = 'k';
        scatter_marker_edgeWidth = 0.5;
        scatter_color_transparency = 0.7; %faceAlpha        

        %% boxplot for each column in Young
        bp1 = boxplot(YoungParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x-shift_value);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %% boxplot for each column in HealthOld
        bp2 = boxplot(HealthyOldParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x);
        set(bp2,'linewidth',box_lineWidth);

        hold on
        %% boxplot for each column in MCIMerge
        bp3 = boxplot(MCIParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x+shift_value);
        set(bp3,'linewidth',box_lineWidth);

        %% Coloring each box
        %findobj first getting the three boxes from right to left
        %first MCI (bp3), then HealthyOld (from bp2), finally Young (from bp1)
        h = findobj(gca,'Tag','Box'); 
        for i = 1:length(h)
            if i<4  %get the MCI box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForMCI,'FaceAlpha',box_color_transparency);
            elseif i<7  %get the HelthyOld box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForHOld,'FaceAlpha',box_color_transparency);
            else %get the Young box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForYoung,'FaceAlpha',box_color_transparency);
            end
        end

        %% Adjusting median
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% add scatter plot and the mean of Young
        num_points = size(YoungParamAllConds,1);
        for i=1:size(YoungParamAllConds,2)
            hold on
            x = i*ones(num_points,1)-shift_value+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, YoungParamAllConds(:,i), scatter_markerSize, ...
                    'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForYoung, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
            %add errorbar
            mean_Young = mean(YoungParamAllConds(:,i));
            sem_Young = std(YoungParamAllConds(:,i))./sqrt(length(YoungParamAllConds(:,i)));
            errorbar(i-shift_value,mean_Young,sem_Young,'k','LineStyle','None', 'LineWidth', 2,  'CapSize', 14); 
            hold on
            %add mean point
            scatter(i-shift_value, mean_Young, 4*scatter_markerSize, 'd',...
                    'filled','MarkerEdgeColor','k', ...
                    'MarkerFaceColor','w', ...
                    'LineWidth',scatter_marker_edgeWidth);
        end

        %% add scatter plot and the mean of HealthyOld
        num_points = size(HealthyOldParamAllConds,1);
        for i=1:size(HealthyOldParamAllConds,2)
            hold on
            x = i*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, HealthyOldParamAllConds(:,i), scatter_markerSize, ...
                    'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForHOld, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
            %add errorbar
            mean_Hold = mean(HealthyOldParamAllConds(:,i));
            sem_Hold = std(HealthyOldParamAllConds(:,i))./sqrt(length(HealthyOldParamAllConds(:,i)));
            errorbar(i,mean_Hold,sem_Hold,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 14);    
            hold on
            %add mean point
            scatter(i, mean_Hold, 4*scatter_markerSize, 'd',...
                    'filled','MarkerEdgeColor','k', ...
                    'MarkerFaceColor','w', ...
                    'LineWidth',scatter_marker_edgeWidth);
        end

        %% add scatter plot and the mean of MCI
        num_points = size(MCIParamAllConds,1);
        for i=1:size(MCIParamAllConds,2)
            hold on
            x = i*ones(num_points,1)+shift_value+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, MCIParamAllConds(:,i), scatter_markerSize, ...
                    'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForMCI, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
            %add errorbar
            mean_MCI = mean(MCIParamAllConds(:,i));
            sem_MCI = std(MCIParamAllConds(:,i))./sqrt(length(MCIParamAllConds(:,i)));
            errorbar(i+shift_value,mean_MCI,sem_MCI,'k','LineStyle','None', 'LineWidth', 3,  'CapSize', 14); 
            hold on
            %add mean point
            scatter(i+shift_value, mean_MCI, 4*scatter_markerSize, 'd',...
                    'filled','MarkerEdgeColor','k', ...
                    'MarkerFaceColor','w', ...
                    'LineWidth',scatter_marker_edgeWidth);
        end

        %% Further post-processing the figure

        %calculate te ylim
        alldata = [YoungParamAllConds;HealthyOldParamAllConds;MCIParamAllConds];
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
            'XTickLabel'  , {'No Change','No Distal Cue', 'No Optical Flow'},...
            'XLim'        , [0.5, 3.5],...
            'YLim'        , lowupYlim,... 
            'LineWidth'   , .5        );
        ylabel(ParamName(ParamIndx));
        allpatches = findall(gca,'type','Patch');
        legend(allpatches(1:3:end), {'Young', 'HealthyOld', 'MCIMerged'}, 'Location','northeast', 'NumColumns',3);

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
        exportgraphics(f,config.ResultFolder+"/Box_"+StoreName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/Box_"+StoreName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end

%%
function plotBoxOfFittedParamMergeCondition(AllYoungParams, AllHealthyOldParams, AllMCIParams, multicomp_tab1, config)
    
    numConds = 3; %3 is the condition number
    ParamName = ["beta", "bG_3", "g_2", "Angular gain g_3", "Anglar bias b", "sigma", "nu"];
    StoreName = ["beta", "bG_3", "g_2", "g_3", "b", "sigma", "nu"];
    for ParamIndx=1:length(ParamName)

        YoungParamAllConds = [];
        HealthyOldParamAllConds = []; %dimension are different, so separate from MCIParamAllConds
        MCIParamAllConds = [];
        
        YoungParamAllCondsMergeinColumn = [];
        HealthyOldParamAllCondsMergeinColumn = [];
        MCIParamAllCondsMergeinColumn = [];

        for TRIAL_FILTER=1:numConds
            %% extract data
            YoungParam = AllYoungParams{TRIAL_FILTER}(:,ParamIndx);
            YoungParamAllConds = [YoungParamAllConds,YoungParam];
            YoungParamAllCondsMergeinColumn = [YoungParamAllCondsMergeinColumn;YoungParam];
            
            HealthyOldParam = AllHealthyOldParams{TRIAL_FILTER}(:,ParamIndx);
            HealthyOldParamAllConds = [HealthyOldParamAllConds,HealthyOldParam];
            HealthyOldParamAllCondsMergeinColumn = [HealthyOldParamAllCondsMergeinColumn;HealthyOldParam];

            MCIParam = AllMCIParams{TRIAL_FILTER}(:,ParamIndx);
            MCIParamAllConds = [MCIParamAllConds,MCIParam];
            MCIParamAllCondsMergeinColumn = [MCIParamAllCondsMergeinColumn;MCIParam];            
        end
    
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
        box_widths_value = 0.4;
        box_color_transparency = 0.5; %faceAlpha
        median_lineWidth = 2;
        median_color = 'k';
        scatter_jitter_value = 0.2;
        scatter_markerSize=10;
        scatter_marker_edgeColor = 'k';
        scatter_marker_edgeWidth = 0.5;
        scatter_color_transparency = 0.7; %faceAlpha        

        hold on
        %% boxplot for each column in Young
        bp1 = boxplot(YoungParamAllCondsMergeinColumn, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 1);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %% boxplot for each column in HealthOld
        bp2 = boxplot(HealthyOldParamAllCondsMergeinColumn, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 2);
        set(bp2,'linewidth',box_lineWidth);

        hold on
        %% boxplot for each column in MCIMerge
        bp3 = boxplot(MCIParamAllCondsMergeinColumn, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 3);
        set(bp3,'linewidth',box_lineWidth);
        

        %% Coloring each box
        %findobj first getting the box for MCI(from bp3) 
        % then for HealthyOld (from bp2) 
        %then getting the box for Young (frm bp1)
        h = findobj(gca,'Tag','Box'); 
        %get the Young box
        patch(get(h(1),'XData'),get(h(1),'YData'),colorForMCI,'FaceAlpha',box_color_transparency);        
        %get the HelthyOld box
        patch(get(h(2),'XData'),get(h(2),'YData'),colorForHOld,'FaceAlpha',box_color_transparency);
        %get the MCI box
        patch(get(h(3),'XData'),get(h(3),'YData'),colorForYoung,'FaceAlpha',box_color_transparency);

        %% Adjusting median
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% add scatter plot and the mean of Young
        num_points = size(YoungParamAllCondsMergeinColumn,1);
        hold on
        x = 1*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, YoungParamAllCondsMergeinColumn, scatter_markerSize, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForYoung, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        hold on
        %add errorbar
        mean_Young = mean(YoungParamAllCondsMergeinColumn);
        sem_Young = std(YoungParamAllCondsMergeinColumn)./sqrt(length(YoungParamAllCondsMergeinColumn));
        errorbar(1,mean_Young,sem_Young,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18); 
        hold on
        %add mean point
        scatter(1, mean_Young, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);    

        %% add scatter plot and the mean of HealthyOld
        num_points = length(HealthyOldParamAllCondsMergeinColumn);
        hold on
        x = 2*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, HealthyOldParamAllCondsMergeinColumn, scatter_markerSize, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForHOld, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        %add errorbar
        mean_Hold = mean(HealthyOldParamAllCondsMergeinColumn);
        sem_Hold = std(HealthyOldParamAllCondsMergeinColumn)./sqrt(length(HealthyOldParamAllCondsMergeinColumn));
        errorbar(2,mean_Hold,sem_Hold,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18);    
        hold on
        %add mean point
        scatter(2, mean_Hold, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);    

        %% add scatter plot and the mean of MCI
        num_points = length(MCIParamAllCondsMergeinColumn);
        hold on
        x = 3*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, MCIParamAllCondsMergeinColumn, scatter_markerSize, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForMCI, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        %add errorbar
        mean_MCI = mean(MCIParamAllCondsMergeinColumn);
        sem_MCI = std(MCIParamAllCondsMergeinColumn)./sqrt(length(MCIParamAllCondsMergeinColumn));
        errorbar(3,mean_MCI,sem_MCI,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18); 
        hold on
        %add mean point
        scatter(3, mean_MCI, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);     

        %% Further post-processing the figure

        %calculate the ylim
        alldata = [YoungParamAllConds;HealthyOldParamAllConds;MCIParamAllConds];
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
            'XTickLabel'  , {'Young','HealthyOld','MCI'},...
            'LineWidth'   , .5        );
        ylabel(ParamName(ParamIndx));

        %extract pvalue for multicomparison of Group effect for showing on the figure
        multicomp_result = multicomp_tab1{ParamIndx};
        % 1       2             3
        % MCI     HealthyOld    Young
        PvalueYoungvsHealthyOld = multicomp_result(3,6); % Young vs. HealthyOld see Two-way anova for details
        PvalueHealthyOldvsMCI = multicomp_result(1,6); % HealthyOld v.s. MCI vs.  see Two-way anova for details
        PvalueYoungvsMCI = multicomp_result(2,6); % Young vs. MCI see Two-way anova for details 
%         str = {['P12 = ',sprintf('%.2g',PvalueYoungvsHealthyOld)],...
%                ['P23 = ',sprintf('%.2g',PvalueHealthyOldvsMCI)],...
%                ['P13 = ',sprintf('%.2g',PvalueYoungvsMCI)]};
%         annotation('textbox',[0.2 0.6 0.3 0.3],'String',str,'FitBoxToText','on');
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
        exportgraphics(f,config.ResultFolder+"/ZMergeCondsBox_"+StoreName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/ZMergeCondsBox_"+StoreName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end

%% Two-way anova on merged MCI
function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnovaOn_Young_HealthyOld_MergeMCI(AllYoung, AllHealthyOld, AllMCI, config)

    %load configurations necessary for the script
    resultfolder = config.ResultFolder;
    
    %create storing folder for trajectory if not exist
    savefoldername = resultfolder+"/TwowayAnova/";
    if ~exist(savefoldername, 'dir')
       mkdir(savefoldername);
    end
    
    param_names = ["beta", "bG3", "g2", "g3", 'b', "sigma", "nu"];
    param_nums = length(param_names);
    
    anova_tab = cell(0);
    multicomp_tab1 = cell(0);
    multicomp_tab2 = cell(0);
    multicomp_tab12 = cell(0);
    %for each parameter, calculate the anova and multiple test
    for param_idx=1:param_nums
        
        param_name = param_names(param_idx);
        
        %processing the data into a long numeric vector 
        [MCIY, MCIGroupNames, MCIConditionNames]=GroupAndRemoveNaN(AllMCI,param_idx,'MCI');
        [HealthyOldY, HealthyOldGroupNames, HealthyOldConditionNames]=GroupAndRemoveNaN(AllHealthyOld,param_idx,'HealthyOld');
        [YoungY, YoungGroupNames, YoungConditionNames]=GroupAndRemoveNaN(AllYoung,param_idx,'Young');
        
        AllY = [MCIY,HealthyOldY,YoungY];
        AllGroupNames = [MCIGroupNames,HealthyOldGroupNames,YoungGroupNames];
        AllConditionNames = [MCIConditionNames,HealthyOldConditionNames,YoungConditionNames];
    
        %Do two-way anova with unbalanced design
        [p,tb1, stats]= anovan(AllY,{AllGroupNames,AllConditionNames},'model','interaction','varnames',{'Groups','Conditions'},'display','on');
        anova_tab{param_idx} = tb1;
    
        %Do multiple comparisons on main effect 1
        result = multcompare(stats,'Dimension',[1],'CType','bonferroni');
        multicomp_tab1{param_idx} = result;
        title("Multiple comparisons with bonferroni correction of : "+param_name);
        saveas(gcf,savefoldername+"MultiCompME1_"+param_name+".png");
        close(gcf);
    
        %Do multiple comparisons on main effect 2
        result = multcompare(stats,'Dimension',[2],'CType','bonferroni');
        multicomp_tab2{param_idx} = result;
        title("Multiple comparisons with bonferroni correction of parameter: "+param_name);
        saveas(gcf,savefoldername+"MultiCompME2_"+param_name+".png");
        close(gcf);
    
        %Do multiple comparisons on main effect 1&2
        result = multcompare(stats,'Dimension',[1,2],'CType','bonferroni');
        multicomp_tab12{param_idx} = result;
        title("Multiple comparisons with bonferroni correction of parameter: "+param_name);
        saveas(gcf,savefoldername+"MultiCompME1ME2_"+param_name+".png");
        close(gcf);    
    end
    
end

