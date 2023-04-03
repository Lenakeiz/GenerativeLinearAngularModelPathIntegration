%% Script to create output for Fig. 4 - parameter comparisons between MCI positive and MCI negative
% Zilong Ji, UCL, 2022 zilong.ji@ucl.ac.uk
% Fits the model on all of the groups and run a two-way Anova
% (group*condition) on MCI positive and MCI negative
% Output: for each parameter fitted by the model output one boxplot with
% the three groups and performance splitted by environmental condition and
% a boxplot with three groups performance averaged across environmental
% conditions

% Preparing the data
GLAMPI_PrepareBaseConfig;

% Preprocessing the data
GLAMPI_PreprocessData;

% Model fitting
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI;

%% Preparing the output
config.ResultFolder     =   pwd + "/Output/Fig4/"+config.ModelName+"/MCIPosvsMCINeg";

if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Generating color scheme for our paper
ColorPattern; 

% Collecting information from output
AllMCIPosParams     =   MCIPos.Results.estimatedParams;
AllMCINegParams     =   MCINeg.Results.estimatedParams;

% TwowayAnova Analysis
[anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnova_MCIPosMCINeg(AllMCIPosParams, AllMCINegParams, config);

% Plot results
BoxPlotOfFittedParam(AllMCIPosParams, AllMCINegParams, anova_tab, config);
BoxPlotOfFittedParamMergeCondition(AllMCIPosParams, AllMCINegParams, multicomp_tab1, config)

% Final cleanup to leave workspace as the end of the Preprocessing stage.
% Remove if you want to take a look at the output data.
clearvars -except config YoungControls HealthyControls MCINeg MCIPos MCIUnk anova_tab multicomp_tab1 multicomp_tab2 multicomp_tab12

%% ---------------------------------------------------------------------
function BoxPlotOfFittedParam(AllMCIPosParams, AllMCINegParams, anova_tab, config)
    
    numConds = 3; % environmental conditions
    ParamName = config.ParamName;
    for ParamIndx=1:length(ParamName)

        MCIPosParamAllConds = []; 
        MCINegParamAllConds = [];
        for TRIAL_FILTER=1:numConds
            %% extract data for each condition
            MCIPosParam = AllMCIPosParams{TRIAL_FILTER}(:,ParamIndx);
            MCIPosParamAllConds = [MCIPosParamAllConds,MCIPosParam];

            MCINegParam = AllMCINegParams{TRIAL_FILTER}(:,ParamIndx);
            MCINegParamAllConds = [MCINegParamAllConds,MCINegParam];
        end

        % remove Nan rows (because of 1) removing participants with short walking length; 2) not enough trials for parameter estimation)
        nonNAN_MCIPosParamAllConds = removeNanRows(MCIPosParamAllConds);
        nonNAN_MCINegParamAllConds = removeNanRows(MCINegParamAllConds);
    
        f = figure('visible','off','Position', [100 100 1000 500]);

        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)     

        colorForMCIPos = config.color_scheme_npg(6,:);
        colorForMCINeg = config.color_scheme_npg(3,:);

        % parameters set for controlling visual output
        whisker_value               =   1.5;
        box_lineWidth               =   0.3;
        box_widths_value            =   0.2;
        box_color_transparency      =   0.5;
        center_x                    =   [1,2,3];    
        %center of box (three conditions)
        shift_value                 =   0.2;        
        %box shift from center
        median_lineWidth            =   2;
        median_color                =   'k';
        scatter_jitter_value        =   0.1;
        scatter_markerSize          =   30;
        scatter_marker_edgeColor    =   'k';
        scatter_marker_edgeWidth    =   0.5;
        scatter_color_transparency  =   0.7;     

        %% Boxplot for each column in MCI positive
        bp1 = boxplot(MCIPosParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ...
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x-shift_value);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        
        %% Boxplot for each column in MCI negative
        bp2 = boxplot(MCINegParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ...
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x+shift_value);
        set(bp2,'linewidth',box_lineWidth);

        %% Boxplot visual changes
        h = findobj(gca,'Tag','Box'); 
        for i = 1:length(h)
            if i<4  %get the MCINeg box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForMCINeg,'FaceAlpha',box_color_transparency);
            else %get the MCIPos box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForMCIPos,'FaceAlpha',box_color_transparency);
            end
        end

        %% Median visual changes
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% Scatter plot for data and mean (MCI positive)
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
            mean_MCIPos = mean(nonNAN_MCIPosParamAllConds(:,i));
            sem_MCIPos = std(nonNAN_MCIPosParamAllConds(:,i))./sqrt(length(nonNAN_MCIPosParamAllConds(:,i)));
            errorbar(i-shift_value,mean_MCIPos,sem_MCIPos,'k','LineStyle','None', 'LineWidth', 2,  'CapSize', 14);    
            hold on
            %add mean point
            scatter(i-shift_value, mean_MCIPos, 4*scatter_markerSize, 'd',...
                    'filled','MarkerEdgeColor','k', ...
                    'MarkerFaceColor','w', ...
                    'LineWidth',scatter_marker_edgeWidth);
        end

        %% Scatter plot for data and mean (MCI negative)
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
            mean_MCI = mean(nonNAN_MCINegParamAllConds(:,i));
            sem_MCI = std(nonNAN_MCINegParamAllConds(:,i))./sqrt(length(nonNAN_MCINegParamAllConds(:,i)));
            errorbar(i+shift_value,mean_MCI,sem_MCI,'k','LineStyle','None', 'LineWidth', 2,  'CapSize', 14); 
            hold on
            %add mean point
            scatter(i+shift_value, mean_MCI, 4*scatter_markerSize, 'd',...
                    'filled','MarkerEdgeColor','k', ...
                    'MarkerFaceColor','w', ...
                    'LineWidth',scatter_marker_edgeWidth);
        end

        %% Figure post-processing
        % calculate the Y limits
        alldata = [MCIPosParamAllConds;MCINegParamAllConds];
        maxdata = max(alldata,[],'all');
        mindata = min(alldata, [], 'all');
        lowupYlim = [mindata-.1*(maxdata-mindata)-eps, maxdata+.1*(maxdata-mindata)+eps]; 

        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.0 .0 .0], ...
            'YColor'      , [.0 .0 .0], ...
            'XTick'       , (1:3),... 
            'XLim'        , [0.5, 3.5],...
            'YLim'        , lowupYlim,...   
            'XTickLabel'  , {'No Change','No Distal Cue', 'No Optical Flow'},...
            'LineWidth'   , 1.0        );
        ylabel(ParamName(ParamIndx));
        allpatches = findall(gca,'type','Patch');
        legend(allpatches(1:3:end), {'MCIPos' 'MCINeg'}, 'Location','northeast', 'NumColumns',2);

        %extract pvalue for group, conditino and interaction
        anova_result = anova_tab{ParamIndx};
        group_pvalue = anova_result{2,7};
        condition_pvalue = anova_result{3,7};
        interaction_pvalue = anova_result{4,7};

        title(strcat(['Group P = ',sprintf('%.2g',group_pvalue)],...
              ['    Condition P = ',sprintf('%.2g',condition_pvalue)],...
              ['    Interaction P = ',sprintf('%.2g',interaction_pvalue)]))

        %% Export figure
        exportgraphics(f,config.ResultFolder+"/Box_"+ParamName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/Box_"+ParamName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end

%% ---------------------------------------------------------------------
function BoxPlotOfFittedParamMergeCondition(AllMCIPosParams, AllMCINegParams, multicomp_tab1, config)
    
    numConds = 3; % environmental conditions
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

        % remove Nan rows (because of 1) removing participants with short walking length; 2) not enough trials for parameter estimation)
        MCIPosParamAllConds = removeNanRows(MCIPosParamAllConds);
        MCINegParamAllConds = removeNanRows(MCINegParamAllConds);

        MCIPosParamMean = mean(MCIPosParamAllConds, 2);
        MCINegParamMean = mean(MCINegParamAllConds, 2);

        f = figure('visible','off','Position', [100 100 500 500]);

        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)     

        colorForMCIPos = config.color_scheme_npg(6,:);
        colorForMCINeg = config.color_scheme_npg(3,:);

        % parameters set for controlling visual output
        whisker_value               =   1.5;
        box_lineWidth               =   0.3;
        box_widths_value            =   0.4;
        box_color_transparency      =   0.5; 
        median_lineWidth            =   2;
        median_color                =   'k';
        scatter_jitter_value        =   0.2;
        scatter_markerSize          =   50;
        scatter_marker_edgeColor    =   'k';
        scatter_marker_edgeWidth    =   0.5;
        scatter_color_transparency  =   0.7;         

        %% Boxplot for each column in MCI positive
        bp1 = boxplot(MCIPosParamMean, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 1);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %% Boxplot for each column in MCI negative
        bp2 = boxplot(MCINegParamMean, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 2);
        set(bp2,'linewidth',box_lineWidth);

        %% Boxplot visual changes
        h = findobj(gca,'Tag','Box'); 
        %get the MCI negative
        patch(get(h(1),'XData'),get(h(1),'YData'),colorForMCINeg,'FaceAlpha',box_color_transparency);
        %get the MCI positive
        patch(get(h(2),'XData'),get(h(2),'YData'),colorForMCIPos,'FaceAlpha',box_color_transparency);

        %% Median visual change
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% Scatter plot for data and mean MCI positive
        num_points = length(MCIPosParamMean);
        hold on
        x = ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5);
        scatter(x, MCIPosParamMean, scatter_markerSize, ...
                'filled', ...
                'o', ... 
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForMCIPos, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        % add errorbar
        mean_MCIPos = mean(MCIPosParamMean);
        sem_MCIPos= std(MCIPosParamMean)./sqrt(length(MCIPosParamMean));
        errorbar(1,mean_MCIPos,sem_MCIPos,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18);    
        hold on
        % add mean point
        scatter(1, mean_MCIPos, 3*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);

        %% Scatter plot for data and mean MCI negative
        num_points = size(MCINegParamMean,1);
        hold on
        x = 2*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5);
        scatter(x, MCINegParamMean, scatter_markerSize, ...
                'filled', ...
                'o', ... 
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
        scatter(2, mean_MCINeg, 3*scatter_markerSize, 'd',...
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
        alldata = [MCIPosParamMean;MCINegParamMean];
        maxdata = max(alldata,[],'all');
        mindata = min(alldata, [], 'all');
        
        if ParamName(ParamIndx)=="g2" | ParamName(ParamIndx)=="g3" | ParamName(ParamIndx)=="nu" 
            lowupYlim = [0, 3];
            yticks = [0,1,2,3];
        elseif ParamName(ParamIndx)=="beta"
            lowupYlim = [-0.1, 0.5];
            yticks = [-0.1, 0, 0.1, 0.3, 0.5];
        elseif ParamName(ParamIndx)=="k"
            lowupYlim = [0, 2];
            yticks = [0,0.5,1,1.5,2.0];  
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
            'XTick'       , (1:2),... 
            'XLim'        , [0.5, 2.5],...
            'YLim'        , lowupYlim,...   
            'YTick'       , yticks,...
            'XTickLabel'  , {'MCIPos','MCINeg'},...
            'LineWidth'   , 1.0        );
        ylabel(ParamName(ParamIndx));

        % Extract pvalues from multiple comparison of group effect. Adding
        % it to the figure
        multicomp_result = multicomp_tab1{ParamIndx};
        Pvalue = multicomp_result(1,6); % MCIPos vs. MCINeg see Two-way anova for details

        title(['P value = ',sprintf('%.2g',Pvalue)])

        %% Add significance bars
        if Pvalue<0.05
            H=adjustablesigstar({[1,2]},[Pvalue]);
        end

        %% export figure
        exportgraphics(f,config.ResultFolder+"/MergeCondsBox_"+ParamName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/MergeCondsBox_"+ParamName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end

