%% Script to create output for Fig. S8 - parameter comparisons between Young, healthy Elderly, MCI unknown, MCI negative, MCI positive
% Andrea Castegnaro, UCL, 2023 uceeaca@ucl.ac.uk
% Fits the model on all of the groups and run a two-way Anova
% (group*condition) on Young Controls, Elderly Controls and MCI Unk/Neg/Pos
% Output: for each parameter fitted by the model output one boxplot with
% the three groups and performance splitted by environmental condition and
% a boxplot with three groups performance averaged across environmental
% conditions

% Preparing the data
GLAMPI_PrepareBaseConfig;

% Preprocessing the data
GLAMPI_PreprocessData;

% Model fitting
config.ModelName        =   "beta_k_g2_g3_m3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "m3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI;

% Generating color scheme for our paper
ColorPattern;

%% Preparing output
config.ResultFolder = pwd + "/Output/FigS9/"+config.ModelName+"/Young_HealthyOld_MCIUnk_MCINeg_MCIPos";
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Collecting information from output
[YoungParams]                    = YoungControls.Results.estimatedParams;
[HealthyParams]                  = HealthyControls.Results.estimatedParams;
[MCIUnkParams]                   = MCIUnk.Results.estimatedParams;
[MCINegParams]                   = MCINeg.Results.estimatedParams;
[MCIPosParams]                   = MCIPos.Results.estimatedParams;

% TwowayAnova Analysis
[anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnova_ModelParams_AllGroups(YoungParams, HealthyParams, MCIUnkParams, MCINegParams, MCIPosParams, config);

BoxPlotOfFittedParamMergeCondition(YoungParams, HealthyParams, MCIUnkParams, MCINegParams, MCIPosParams,multicomp_tab1, config)

% Final cleanup to leave workspace as the end of the Preprocessing stage.
% Remove if you want to take a look at the output data.
clearvars -except config YoungControls HealthyControls MCINeg MCIPos MCIUnk anova_tab multicomp_tab1 multicomp_tab2 multicomp_tab12

%% ---------------------------------------------------------------------
function BoxPlotOfFittedParamMergeCondition(AllYoungParams, AllHealthyOldParams, AllMCIUnkParams, AllMCINegParams, AllMCIPosParams, multicomp_tab1, config)
    
    numConds = 3; % environmental conditions
    ParamName = config.ParamName;
    for ParamIndx=1:length(ParamName)

        YoungParamAllConds = [];
        HealthyOldParamAllConds = []; 
        MCIUnkParamAllConds = [];
        MCINegParamAllConds = [];
        MCIPosParamAllConds = [];
        
        for TRIAL_FILTER=1:numConds
            %% extract data
            YoungParam          = AllYoungParams{TRIAL_FILTER}(:,ParamIndx);
            YoungParamAllConds  = [YoungParamAllConds,YoungParam];
            
            HealthyOldParam     = AllHealthyOldParams{TRIAL_FILTER}(:,ParamIndx);
            HealthyOldParamAllConds = [HealthyOldParamAllConds,HealthyOldParam];

            MCIUnkParams            = AllMCIUnkParams{TRIAL_FILTER}(:,ParamIndx);
            MCIUnkParamAllConds    = [MCIUnkParamAllConds,MCIUnkParams];

            MCINegParams            = AllMCINegParams{TRIAL_FILTER}(:,ParamIndx);
            MCINegParamAllConds    = [MCINegParamAllConds,MCINegParams]; 

            MCIPosParams            = AllMCIPosParams{TRIAL_FILTER}(:,ParamIndx);
            MCIPosParamAllConds    = [MCIPosParamAllConds,MCIPosParams]; 
        end

        % remove Nan rows (because of 1) removing participants with short walking length; 2) not enough trials for parameter estimation)
        YoungParamAllConds      = removeNanRows(YoungParamAllConds);
        HealthyOldParamAllConds = removeNanRows(HealthyOldParamAllConds);
        MCIUnkParamAllConds     = removeNanRows(MCIUnkParamAllConds);
        MCINegParamAllConds     = removeNanRows(MCINegParamAllConds);
        MCIPosParamAllConds     = removeNanRows(MCIPosParamAllConds);
        
        YoungParamMean      = mean(YoungParamAllConds, 2);
        HealthyOldParamMean = mean(HealthyOldParamAllConds, 2);
        MCIUnkParamMean        = mean(MCIUnkParamAllConds, 2);
        MCINegParamMean        = mean(MCINegParamAllConds, 2);
        MCIPosParamMean        = mean(MCIPosParamAllConds, 2);
    
        f = figure('visible','on','Position', [100 100 360 360]);

        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)     

        colorForYoung = config.color_scheme_npg(1,:);        
        colorForElderly = config.color_scheme_npg(5,:);
        colorForMCIUnk = config.color_scheme_npg(2,:);
        colorForMCINeg = config.color_scheme_npg(4,:);
        colorForMCIPos = config.color_scheme_npg(6,:);

        % parameters set for controlling visual output
        whisker_value               =   1.5;
        box_lineWidth               =   0.3;
        box_widths_value            =   0.4;
        box_color_transparency      =   0.5; %faceAlpha
        median_lineWidth            =   2;
        median_color                =   'k';
        scatter_jitter_value        =   0.2;
        scatter_markerSize          =   50;
        scatter_marker_edgeColor    =   'k';
        scatter_marker_edgeWidth    =   0.5;
        scatter_color_transparency  =   0.7; %faceAlpha        

        hold on
        %% Boxplot for each column in Young
        bp1 = boxplot(YoungParamMean, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 1);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %% Boxplot for each column in Healthy Elderly
        bp2 = boxplot(HealthyOldParamMean, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 2);
        set(bp2,'linewidth',box_lineWidth);

        hold on
        %% Boxplot for each column in MCI negative
        bp3 = boxplot(MCIUnkParamMean, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 3);
        set(bp3,'linewidth',box_lineWidth);

                hold on
        %% Boxplot for each column in MCI negative
        bp4 = boxplot(MCINegParamMean, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 4);
        set(bp4,'linewidth',box_lineWidth);

                hold on
        %% Boxplot for each column in MCI negative
        bp5 = boxplot(MCIPosParamMean, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 5);
        set(bp5,'linewidth',box_lineWidth);

        %% Boxplot visual changes
        % matlab has a lifo system for these
        h = findobj(gca,'Tag','Box'); 
        %get the Young box
        patch(get(h(5),'XData'),get(h(5),'YData'),colorForYoung,'FaceAlpha',box_color_transparency);        
        %get the HelthyOld box
        patch(get(h(4),'XData'),get(h(4),'YData'),colorForElderly,'FaceAlpha',box_color_transparency);
        %get the MCI Unk box
        patch(get(h(3),'XData'),get(h(3),'YData'),colorForMCIUnk,'FaceAlpha',box_color_transparency);
        %get the MCI Unk box
        patch(get(h(2),'XData'),get(h(2),'YData'),colorForMCINeg,'FaceAlpha',box_color_transparency);
        %get the MCI Unk box
        patch(get(h(1),'XData'),get(h(1),'YData'),colorForMCIPos,'FaceAlpha',box_color_transparency);

        %% Median visual change
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% Scatter plot for data and mean (Young)
        num_points = size(YoungParamMean,1);
        hold on
        x = 1*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, YoungParamMean, scatter_markerSize, ...
                'filled', ...
                'o', ... 
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForYoung, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        hold on
        %add errorbar
        mean_Young = mean(YoungParamMean);
        sem_Young = std(YoungParamMean)./sqrt(length(YoungParamMean));
        errorbar(1,mean_Young,sem_Young,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18); 
        hold on
        %add mean point
        scatter(1, mean_Young, 3*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);    

        %% Scatter plot for data and mean (Healthy Elderly)
        num_points = length(HealthyOldParamMean);
        hold on
        x = 2*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, HealthyOldParamMean, scatter_markerSize, ...
                'filled', ...
                'o', ...
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForElderly, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        %add errorbar
        mean_Hold = mean(HealthyOldParamMean);
        sem_Hold = std(HealthyOldParamMean)./sqrt(length(HealthyOldParamMean));
        errorbar(2,mean_Hold,sem_Hold,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18);    
        hold on
        %add mean point
        scatter(2, mean_Hold, 3*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);    

        %% Scatter plot for data and mean (MCI unknown)
        num_points = length(MCIUnkParamMean(:,:));
        hold on
        x = 3*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, MCIUnkParamMean(:,:), scatter_markerSize, ...
                'filled', ...
                'square', ...
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForMCIUnk, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth);

        %add mean + errorbar
        mean_MCI = mean(MCIUnkParamMean);
        sem_MCI = std(MCIUnkParamMean)./sqrt(length(MCIUnkParamMean));
        errorbar(3,mean_MCI,sem_MCI,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18); 
        hold on
        %add mean point
        scatter(3, mean_MCI, 3*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);     

        %% Scatter plot for data and mean (MCI negative)
        num_points = length(MCINegParamMean(:,:));
        hold on
        x = 4*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, MCINegParamMean(:,:), scatter_markerSize, ...
                'filled', ...
                'v', ...
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForMCINeg, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth);

        %add mean + errorbar
        mean_MCI = mean(MCINegParamMean);
        sem_MCI = std(MCINegParamMean)./sqrt(length(MCINegParamMean));
        errorbar(4,mean_MCI,sem_MCI,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18); 
        hold on
        %add mean point
        scatter(4, mean_MCI, 3*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);     

        %% Scatter plot for data and mean (MCI positive)
        num_points = length(MCIPosParamMean(:,:));
        hold on
        x = 5*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, MCIPosParamMean(:,:), scatter_markerSize, ...
                'filled', ...
                '^', ...
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForMCIPos, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth);

        %add mean + errorbar
        mean_MCI = mean(MCIPosParamMean);
        sem_MCI = std(MCIPosParamMean)./sqrt(length(MCIPosParamMean));
        errorbar(5,mean_MCI,sem_MCI,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18); 
        hold on
        %add mean point
        scatter(5, mean_MCI, 3*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);     
        
        % add horizontal line for parameter reference
        if ParamName(ParamIndx)=="beta"
            yline(0,Color='r',LineStyle='--',LineWidth=2);
        elseif ParamName(ParamIndx)=="k"
            yline(1,Color='r',LineStyle='--',LineWidth=2);
        elseif ParamName(ParamIndx)=="g2"
            yline(1,Color='r',LineStyle='--',LineWidth=2);
        elseif ParamName(ParamIndx)=="g3"
            yline(1,Color='r',LineStyle='--',LineWidth=2);
        end

        %% Figure post-processing
        % calculate the Y limits
        alldata = [YoungParamMean;HealthyOldParamMean;MCIUnkParamMean;MCINegParamMean;MCIPosParamMean];
        maxdata = max(alldata,[],'all');
        mindata = min(alldata, [], 'all');

        if ParamName(ParamIndx)=="g2" | ParamName(ParamIndx)=="g3" | ParamName(ParamIndx)=="nu" 
            lowupYlim = [0, 3];
            yticks = [0,1,2,3];
        elseif ParamName(ParamIndx)=="k"
            lowupYlim = [0, 2];
            yticks = [0,0.5,1,1.5,2.0];            
        elseif ParamName(ParamIndx)=="beta"
            lowupYlim = [-0.1, 0.5];
            yticks = [-0.1, 0, 0.1, 0.3, 0.5];
        elseif ParamName(ParamIndx)=="sigma"
            lowupYlim = [0, 1.5];    
            yticks = [0, 0.5, 1.0, 1.5];
        end

        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.0 .0 .0], ...
            'YColor'      , [.0 .0 .0], ...
            'XTick'       , (1:5),... 
            'XLim'        , [0.5, 5.5],...
            'YLim'        , lowupYlim,...
            'YTick'       , yticks,...
            'XTickLabel'  , {'Young','Elderly','MCI Unk','MCI Neg','MCI Pos'},...
            'LineWidth'   , 1.0        );

        ylabel(ParamName(ParamIndx));

        % Extract pvalues from multiple comparison of group effect. Adding
        % it to the figure
        multicomp_result = multicomp_tab1{ParamIndx};
        % 1       2             3        4         5
        % Young   HealthyOld    MCI Unk  MCI Neg   MCI Pos
        PvalueYoungvsHealthyOld = multicomp_result(1,6); % Young vs. HealthyOld
        PvalueHealthyOldvsMCIUnk = multicomp_result(5,6); % HealthyOld v.s. MCI unk
        PvalueHealthyOldvsMCINeg = multicomp_result(6,6); % HealthyOld v.s. MCI neg
        PvalueHealthyOldvsMCIPos = multicomp_result(7,6); % HealthyOld v.s. MCI pos
        PvalueMCINegvsMCIPos     = multicomp_result(10,6);
        PvalueMCIUnkvsMCIPos     = multicomp_result(9,6); % MCI Unk vs MCI pos
        PvalueMCIUnkvsMCINeg     = multicomp_result(8,6); % MCI Unk vs MCI neg

        
        %% Add significance bars
        AllP = [PvalueYoungvsHealthyOld,PvalueHealthyOldvsMCIUnk,PvalueHealthyOldvsMCINeg,PvalueHealthyOldvsMCIPos,PvalueMCINegvsMCIPos,PvalueMCIUnkvsMCIPos,PvalueMCIUnkvsMCINeg];
        Xval = [[1,2];[2,3];[2,4];[2,5];[4,5];[3 5];[3 4]];
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
            H=adjustablesigstar(XXX,AllP_select);
        end

        %% export figure
        exportgraphics(f,config.ResultFolder+"/MergeCondsBox_"+ParamName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/MergeCondsBox_"+ParamName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end
