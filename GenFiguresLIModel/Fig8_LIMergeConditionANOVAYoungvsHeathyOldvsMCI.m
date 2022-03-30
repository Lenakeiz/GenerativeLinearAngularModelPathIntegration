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

config.UseGlobalSearch = true;

resultfolder = savefolder+"PaperFigs/Fig8F_GammaFull";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% Model fitting using LeakyIntegration Model
config.ModelName = "GammaFull";
config.NumParams = 5;

%% merge three conditions
config.TrialFilter = 0;

%% Model fitting for YoungControl data
% calculating tracking path and transoform data
config.Speed.tresholdForBadParticipantL1Recontruction = 1.55;   % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
YoungControls   = CalculateTrackingPath(YoungControls, config);
%transform data
YoungControls = TransformPaths(YoungControls);
ResultsYoung = PerformGroupFit(YoungControls, config);
YoungParams = ResultsYoung.estimatedParams;      

%% Model fitting for HealthyOld data
config.Speed.tresholdForBadParticipantL1Recontruction = 2.0; 
HealthyControls   = CalculateTrackingPath(HealthyControls, config);
%transform data
HealthyOld = TransformPaths(HealthyControls);
ResultsHealthyOld = PerformGroupFit(HealthyOld, config);
HealthyOldParams = ResultsHealthyOld.estimatedParams;  

%% Model fitting for MCIPos
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCIPos   = CalculateTrackingPath(MCIPos, config);
%transform data
MCIPos = TransformPaths(MCIPos);
ResultsMCIPos = PerformGroupFit(MCIPos, config);
MCIPosParams = ResultsMCIPos.estimatedParams;  

%% Model fitting for MCINeg
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
MCINeg   = CalculateTrackingPath(MCINeg, config);
%transform data
MCINeg = TransformPaths(MCINeg);
ResultsMCINeg = PerformGroupFit(MCINeg, config);
MCINegParams = ResultsMCINeg.estimatedParams; 

%% Model fitting for MCIUnk
config.Speed.tresholdForBadParticipantL1Recontruction = 0.0; 
Unknown   = CalculateTrackingPath(Unknown, config);
%transform data
MCIUnk = TransformPaths(Unknown);
ResultsMCIUnk = PerformGroupFit(MCIUnk, config);
MCIUnkParams = ResultsMCIUnk.estimatedParams; 

%% Setting colors for using in plots
ColorPattern; 

%% merge MCI together
MCIParams = [MCIPosParams;MCINegParams;MCIUnkParams];

%% OnewayAnova
[~,multicomp_tab] = OnewayAnovaOn_Young_HealthyOld_MergeMCI(YoungParams, HealthyOldParams, MCIParams, config);

%% BarScatter Plot between Young and HealthyOld for all Fitted Params
plotBoxOfFittedParam(YoungParams, HealthyOldParams, MCIParams, multicomp_tab, config);

%%
function plotBoxOfFittedParam(YoungParams, HealthyOldParams, MCIParams, multicomp_tab, config)
    
    ParamName = ["gamma", "bG_3", "g_2", "Angular gain g_3", "Anglar bias b", "sigma", "nu"];
    StoreName = ["gamma", "bG_3", "g_2", "g_3", "b", "sigma", "nu"];
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
        alldata = [MCIParams(:,ParamIndx); HealthyOldParams(:,ParamIndx); YoungParams(:,ParamIndx)];
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
            'XTick'       , (1:3),... 
            'XLim'        , [0.5, 3.5],...
            'YLim'        , lowupYlim,...
            'XTickLabel'  , {'Young','HealthyOld','MCI'},...
            'LineWidth'   , .5        );
            %'Ytick'       , [0,0.5,1.0,1.5],...
            %'YLim'        , [0, 1.5],...   

        ylabel(ParamName(ParamIndx));

        %extract pvalue for multicomparison of Group effect for showing on the figure
        multicomp_result = multicomp_tab{ParamIndx};
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

%% One-way anova on merged MCI
function [anova_tab,multicomp_tab] = OnewayAnovaOn_Young_HealthyOld_MergeMCI(YoungParams, HealthyOldParams, MCIParams, config)

    %load configurations necessary for the script
    resultfolder = config.ResultFolder;
    
    %create storing folder for trajectory if not exist
    savefoldername = resultfolder+"/OnewayAnova/";
    if ~exist(savefoldername, 'dir')
       mkdir(savefoldername);
    end
    
    param_names = ["gamma", "bG3", "g2", "g3", 'b', "sigma", "nu"];
    param_nums = length(param_names);
    
    anova_tab = cell(0);
    multicomp_tab = cell(0);

    %for each parameter, calculate the anova and multiple test
    for param_idx=1:param_nums
        
        param_name = param_names(param_idx);
        
        %processing the data into a long numeric vector 
        Params = [MCIParams(:,param_idx)',HealthyOldParams(:,param_idx)',YoungParams(:,param_idx)'];
        GroupNames = [string(repmat({'MCI'},1,size(MCIParams,1))),...
                      string(repmat({'HealthyOld'},1,size(HealthyOldParams,1))),...
                      string(repmat({'Young'},1,size(YoungParams,1)))];
    
        %Do one-way anova with unbalanced design
        [p,tb1, stats]= anova1(Params, GroupNames, 'display','on');
        anova_tab{param_idx} = tb1;
    
        %Do multiple comparisons on main effect 1
        result = multcompare(stats,'CType','bonferroni');
        multicomp_tab{param_idx} = result;
        title("Multiple comparisons with bonferroni correction of : "+param_name);
        saveas(gcf,savefoldername+"MultiComp"+param_name+".png");
        close(gcf);
    
    end
    
end

