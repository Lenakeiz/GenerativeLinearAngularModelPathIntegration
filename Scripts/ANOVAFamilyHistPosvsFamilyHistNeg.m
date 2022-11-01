%% Cleaning variables and set intial seed for code reproducibility
clearvars; close all; clc;
rng('default'); 

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('AllDataErrorsPrevent.mat');

%% setting the configuration
config.Speed.alpha                                      = 0.9;    % Paramanter for running speed calculation
config.Speed.timeOffsetAfterFlagReach                   = 1.5;    % Time to track after flag reached in seconds 
config.Speed.smoothWindow                               = 10;     % tracking rate should be 10Hz so 4 secs window is 40 datapoints
config.Speed.velocityCutoff                             = 0.2;    % velocity cutoff to select only the walking part of the reconstructed velocity
config.Speed.timeOffsetForDetectedTemporalWindow        = 0.4;    % time in seconds that will push earlier/ the detected rising edge
config.UseGlobalSearch                                  = true;
config.TrackedInboundAngularDeltaT                      = 1;
config.includeStand                                     = false;
config.useweber                                         = false;  % only true when use weber law in simple generative models
config.Speed.tresholdForBadParticipantL1Recontruction   = 0.0;    % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
config.useOoBtrials                                     = true;
config.useTrialFilter                                   = true;

%% Model fitting
%Model related parameters
% config.ModelName        = "beta_g3_sigma_nu";
% config.ParamName        = ["beta", "g3", "sigma", "nu"];

% config.ModelName        = "beta_g2_sigma_nu";
% config.ParamName        = ["beta", "g2", "sigma", "nu"];

config.ModelName        = "beta_g2_g3_sigma_nu";
config.ParamName        = ["beta", "g2", "g3", "sigma", "nu"];

% config.ModelName        = "beta_g2_g3_k3_sigma_nu";
% config.ParamName        = ["beta", "g2", "g3", "k3", "sigma", "nu"];

config.NumParams        = length(config.ParamName);

resultfolder = pwd + "/Output/ModelFigures/FamilyHistPos_FamilyHistNeg_"+config.ModelName;
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% Model fitting for Pos data
FamilyHistPos   = TransformPaths(FamilyHistPos);%transform data
FamilyHistPos   = CalculateTrackingPath(FamilyHistPos, config);
ManuallyScoringFamilyHistPos;
%%
FamilyHistPos.Results = getResultsAllConditions(FamilyHistPos, config);

%% Model fitting for Neg data
FamilyHistNeg   = TransformPaths(FamilyHistNeg);%transform data
FamilyHistNeg   = CalculateTrackingPath(FamilyHistNeg, config);
ManuallyScoringFamilyHistNeg;
%%
FamilyHistNeg.Results = getResultsAllConditions(FamilyHistNeg, config);

%%
AllFamilyHistPosParams     =   FamilyHistPos.Results.estimatedParams;
AllFamilyHistNegParams     =   FamilyHistNeg.Results.estimatedParams;

%% Setting colors for using in plots
ColorPattern; 

%% TwowayAnova and BarScatter Plot
[anova_tab,multicomp_tab1,~, ~] = TwowayAnova_CocoData(AllFamilyHistPosParams, AllFamilyHistNegParams, config);
BoxPlotOfFittedParam(AllFamilyHistPosParams, AllFamilyHistNegParams, anova_tab, config);
BoxPlotOfFittedParamMergeCondition(AllFamilyHistPosParams, AllFamilyHistNegParams, multicomp_tab1, config)

%% ThreewayAnova
config.Gender_Pos = FamilyHistPos.Gender;
config.Gender_Neg = FamilyHistNeg.Gender;
[anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab3, multicomp_tab12] = ThreewayAnova_CocoData(AllFamilyHistPosParams, AllFamilyHistNegParams, config);

%%
function BoxPlotOfFittedParam(AllPosParams, AllNegParams, anova_tab, config)
    
    numConds = 3; %3 is the condition number
    ParamName = config.ParamName;
    for ParamIndx=1:length(ParamName)

        PosParamAllConds = [];
        NegParamAllConds = []; %dimension are different, so separate from YoungParams
        
        for TRIAL_FILTER=1:numConds
            %% extract data
            PosParam            =   AllPosParams{TRIAL_FILTER}(:,ParamIndx);
            PosParamAllConds    =   [PosParamAllConds,PosParam];

            NegParam            =   AllNegParams{TRIAL_FILTER}(:,ParamIndx);
            NegParamAllConds    =   [NegParamAllConds,NegParam];         
        end
        
        %remove Nan rows (nan coz of 1, removing participants with short walking length; 2, not enough trials for parameter estimation)
        PosParamAllConds        =   removeNanRows(PosParamAllConds);
        NegParamAllConds        =   removeNanRows(NegParamAllConds);
    
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
        colorForPos = config.color_scheme_npg(3,:);
        colorForNeg = config.color_scheme_npg(5,:);
        
        %set params
        whisker_value = 1.5;
        box_lineWidth = 0.3;
        box_widths_value = 0.2;
        box_color_transparency = 0.5; %faceAlpha
        %center of box (three conditions)
        center_x = [1,2,3];
        shift_value = 0.15; %box shift from center
        median_lineWidth = 2;
        median_color = 'k';
        scatter_jitter_value = 0.1;
        scatter_markerSize=10;
        scatter_marker_edgeColor = 'k';
        scatter_marker_edgeWidth = 0.5;
        scatter_color_transparency = 0.7; %faceAlpha        

        %% boxplot for each column in Pos
        bp1 = boxplot(PosParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x-shift_value);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %% boxplot for each column in Neg
        bp2 = boxplot(NegParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x+shift_value);
        set(bp2,'linewidth',box_lineWidth);

        hold on

        %% Coloring each box
        %findobj first getting the three boxes for Neg (from bp2) from right to left
        %then getting the three boxes for Pos (frm bp1) from right to left
        h = findobj(gca,'Tag','Box'); 
        for i = 1:length(h)
            if i<4  %get the Neg box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForNeg,'FaceAlpha',box_color_transparency);
            else %get the Pos box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForPos,'FaceAlpha',box_color_transparency);
            end
        end

        %% Adjusting median
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% add scatter plot and the mean of Young
        num_points = size(PosParamAllConds,1);
        for i=1:size(PosParamAllConds,2)
            hold on
            x = i*ones(num_points,1)-shift_value+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, PosParamAllConds(:,i), scatter_markerSize, ...
                    'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForPos, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
            %add errorbar
            mean_Pos = mean(PosParamAllConds(:,i));
            sem_Pos = std(PosParamAllConds(:,i))./sqrt(length(PosParamAllConds(:,i)));
            errorbar(i-shift_value,mean_Pos,sem_Pos,'k','LineStyle','None', 'LineWidth', 2,  'CapSize', 14); 
            hold on
            %add mean point
            scatter(i-shift_value, mean_Pos, 4*scatter_markerSize, 'd',...
                    'filled','MarkerEdgeColor','k', ...
                    'MarkerFaceColor','w', ...
                    'LineWidth',scatter_marker_edgeWidth);
        end

        %% add scatter plot and the mean of Neg
        num_points = size(NegParamAllConds,1);
        for i=1:size(NegParamAllConds,2)
            hold on
            x = i*ones(num_points,1)+shift_value+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, NegParamAllConds(:,i), scatter_markerSize, ...
                    'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForNeg, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
            %add errorbar
            mean_Neg = mean(NegParamAllConds(:,i));
            sem_Neg = std(NegParamAllConds(:,i))./sqrt(length(NegParamAllConds(:,i)));
            errorbar(i+shift_value,mean_Neg,sem_Neg,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 14);    
            hold on
            %add mean point
            scatter(i+shift_value, mean_Neg, 4*scatter_markerSize, 'd',...
                    'filled','MarkerEdgeColor','k', ...
                    'MarkerFaceColor','w', ...
                    'LineWidth',scatter_marker_edgeWidth);
        end

        %% Further post-processing the figure

        %calculate te ylim
        alldata = [PosParamAllConds;NegParamAllConds];
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
        legend(allpatches(1:3:end), {'Pos', 'Neg'}, 'Location','northeast', 'NumColumns',3);

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
function BoxPlotOfFittedParamMergeCondition(AllPosParams, AllNegParams, multicomp_tab1, config)
    
    numConds = 3; %3 is the condition number
    ParamName = config.ParamName;
    for ParamIndx=1:length(ParamName)

        PosParamAllConds = []; 
        NegParamAllConds = [];

        for TRIAL_FILTER=1:numConds
            %% extract data
            PosParam = AllPosParams{TRIAL_FILTER}(:,ParamIndx);
            PosParamAllConds = [PosParamAllConds,PosParam];

            NegParam = AllNegParams{TRIAL_FILTER}(:,ParamIndx);
            NegParamAllConds = [NegParamAllConds,NegParam];
        end

        %remove Nan rows (nan coz of 1, removing participants with short walking length; 2, not enough trials for parameter estimation)
        PosParamAllConds = removeNanRows(PosParamAllConds);
        NegParamAllConds = removeNanRows(NegParamAllConds);

        PosParamMean = mean(PosParamAllConds, 2);
        NegParamMean = mean(NegParamAllConds, 2);
    
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
        colorForPos = config.color_scheme_npg(6,:);
        colorForNeg = config.color_scheme_npg(3,:);

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

        %% boxplot for each column in Pos
        bp1 = boxplot(PosParamMean, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 1);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %boxplot for each column in Neg
        bp2 = boxplot(NegParamMean, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 2);
        set(bp2,'linewidth',box_lineWidth);

        %% Coloring each box
        %findobj first getting the box for Neg (from bp2) 
        %then getting the box for Pos(frm bp1)
        h = findobj(gca,'Tag','Box'); 
        %get the Pos box
        patch(get(h(1),'XData'),get(h(1),'YData'),colorForNeg,'FaceAlpha',box_color_transparency);
        %get the Neg box
        patch(get(h(2),'XData'),get(h(2),'YData'),colorForPos,'FaceAlpha',box_color_transparency);

        %% Adjusting median
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% add scatter plot and the mean of Pos
        num_points = length(PosParamMean);
        hold on
        x = ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, PosParamMean, scatter_markerSize, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForPos, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        %add errorbar
        mean_Pos = mean(PosParamMean);
        sem_Pos= std(PosParamMean)./sqrt(length(PosParamMean));
        errorbar(1,mean_Pos,sem_Pos,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18);    
        hold on
        %add mean point
        scatter(1, mean_Pos, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);

        %% add scatter plot and the mean of Neg
        num_points = size(NegParamMean,1);
        hold on
        x = 2*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, NegParamMean, scatter_markerSize, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForNeg, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        %add errorbar
        mean_Neg = mean(NegParamMean);
        sem_Neg = std(NegParamMean)./sqrt(length(NegParamMean));
        errorbar(2,mean_Neg,sem_Neg,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18); 
        hold on
        %add mean point
        scatter(2, mean_Neg, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);        

        %% Further post-processing the figure
        %calculate the ylim
        alldata = [PosParamMean;NegParamMean];
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
            'XTickLabel'  , {'Pos','Neg'},...
            'LineWidth'   , .5        );
        ylabel(ParamName(ParamIndx));

        %extract pvalue for multicomparison of Group effect for showing on the figure
        multicomp_result = multicomp_tab1{ParamIndx};
        Pvalue = multicomp_result(1,6); % Pos vs. Neg see Two-way anova for details
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