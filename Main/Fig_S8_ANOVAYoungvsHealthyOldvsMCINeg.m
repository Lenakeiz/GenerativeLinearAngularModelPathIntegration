%% Script to create output for Fig. S8 - parameter comparisons between Young, healthy Elderly, MCI negative
% Zilong Ji, UCL, 2022 zilong.ji@ucl.ac.uk
% Fits the model on all of the groups and run a two-way Anova
% (group*condition) on Young Controls, Elderly Controls and MCI negative
% Output: for each parameter fitted by the model output one boxplot with
% the three groups and performance splitted by environmental condition and
% a boxplot with three groups performance averaged across environmental
% conditions

% Preparing the data
VAM_PrepareBaseConfig;

% Preprocessing the data
VAM_PreprocessData;

% Model fitting
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

VAM;

% Preparing output
config.ResultFolder     =   pwd + "/Output/FigS8/"+config.ModelName+"/Young_HealthyOld_MCINeg";
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Generating color scheme for our paper
ColorPattern; 

% Collecting information from output
AllYoungParams      =   YoungControls.Results.estimatedParams;
AllHealthyOldParams =   HealthyControls.Results.estimatedParams;
AllMCINegParams     =   MCINeg.Results.estimatedParams;

% TwowayAnova Analysis
[anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnova_YoungHealthyOldMCICombined(AllYoungParams, AllHealthyOldParams, AllMCINegParams, config);

% Plot results
BoxPlotOfFittedParam(AllYoungParams, AllHealthyOldParams, AllMCINegParams, anova_tab, config);
BoxPlotOfFittedParamMergeCondition(AllYoungParams, AllHealthyOldParams, AllMCINegParams, multicomp_tab1, config)

%%
function BoxPlotOfFittedParam(AllYoungParams, AllHealthyOldParams, AllMCIParams, anova_tab, config)
    
    numConds = 3; % environmental conditions
    ParamName = config.ParamName;
    for ParamIndx=1:length(ParamName)

        YoungParamAllConds = [];
        HealthyOldParamAllConds = []; 
        MCIParamAllConds = [];
        
        for TRIAL_FILTER=1:numConds
            %% extract data
            YoungParam          = AllYoungParams{TRIAL_FILTER}(:,ParamIndx);
            YoungParamAllConds  = [YoungParamAllConds,YoungParam];

            HealthyOldParam     = AllHealthyOldParams{TRIAL_FILTER}(:,ParamIndx);
            HealthyOldParamAllConds = [HealthyOldParamAllConds,HealthyOldParam];

            MCIParam            = AllMCIParams{TRIAL_FILTER}(:,ParamIndx);
            MCIParamAllConds    = [MCIParamAllConds,MCIParam];            
        end
        
        % remove Nan rows (because of 1) removing participants with short walking length; 2) not enough trials for parameter estimation)
        YoungParamAllConds      = removeNanRows(YoungParamAllConds);
        HealthyOldParamAllConds = removeNanRows(HealthyOldParamAllConds);
        MCIParamAllConds        = removeNanRows(MCIParamAllConds);

        f = figure('visible','off','Position', [100 100 1000 500]);

        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)     

        colorForYoung = config.color_scheme_npg(3,:);
        colorForHOld = config.color_scheme_npg(5,:);
        colorForMCI = config.color_scheme_npg(2,:);
        
        % parameters set for controlling visual output
        whisker_value = 1.5;
        box_lineWidth = 0.3;
        box_widths_value = 0.2;
        box_color_transparency = 0.5;
        %center of box (three conditions)
        center_x = [1,2,3];
        shift_value = 0.25; 
        %box shift from center
        median_lineWidth = 2;
        median_color = 'k';
        scatter_jitter_value = 0.1;
        scatter_markerSize=30;
        scatter_marker_edgeColor = 'k';
        scatter_marker_edgeWidth = 0.5;
        scatter_color_transparency = 0.7;         

        %% Boxplot for each column in Young
        bp1 = boxplot(YoungParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x-shift_value);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %% Boxplot for each column in Health Elderly
        bp2 = boxplot(HealthyOldParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x);
        set(bp2,'linewidth',box_lineWidth);

        hold on
        %% Boxplot for each column in MCI Negative
        bp3 = boxplot(MCIParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x+shift_value);
        set(bp3,'linewidth',box_lineWidth);

       %% Boxplot visual changes
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

        %% Median visual changes
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% Scatter plot for data and mean (Young)
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

        %% Scatter plot for data and mean (Healthy Elderly)
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

        %% Scatter plot for data and mean (MCI negative)
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

        %% Figure post-processing
        % calculate the Y limits
        alldata = [YoungParamAllConds;HealthyOldParamAllConds;MCIParamAllConds];
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
            'XTickLabel'  , {'No Change','No Distal Cue', 'No Optical Flow'},...
            'XLim'        , [0.5, 3.5],...
            'YLim'        , lowupYlim,... 
            'LineWidth'   , 1.0        );
        ylabel(ParamName(ParamIndx));
        allpatches = findall(gca,'type','Patch');
        legend(allpatches(1:3:end), {'Young', 'HealthyOld', 'MCINeg'}, 'Location','northeast', 'NumColumns',3);

        %extract pvalue for group, environmental condition and interaction
        anova_result = anova_tab{ParamIndx};
        group_pvalue = anova_result{2,7};
        condition_pvalue = anova_result{3,7};
        interaction_pvalue = anova_result{4,7};

        title(strcat(['Group P = ',sprintf('%.2g',group_pvalue)],...
              ['    Condition P = ',sprintf('%.2g',condition_pvalue)],...
              ['    Interaction P = ',sprintf('%.2g',interaction_pvalue)]))

        %% save figure
        exportgraphics(f,config.ResultFolder+"/Box_"+ParamName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/Box_"+ParamName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end

%% Box Plot Of Fitted Parameters by avearaging three paramters from three conditions into one mean value
function BoxPlotOfFittedParamMergeCondition(AllYoungParams, AllHealthyOldParams, AllMCIParams,multicomp_tab1, config)
    
    numConds = 3; % environmental conditions
    ParamName = config.ParamName;
    for ParamIndx=1:length(ParamName)

        YoungParamAllConds = [];
        HealthyOldParamAllConds = []; 
        MCIParamAllConds = [];
        for TRIAL_FILTER=1:numConds
            %% extract data
            YoungParam          = AllYoungParams{TRIAL_FILTER}(:,ParamIndx);
            YoungParamAllConds  = [YoungParamAllConds,YoungParam];
            
            HealthyOldParam     = AllHealthyOldParams{TRIAL_FILTER}(:,ParamIndx);
            HealthyOldParamAllConds = [HealthyOldParamAllConds,HealthyOldParam];

            MCIParam            = AllMCIParams{TRIAL_FILTER}(:,ParamIndx);
            MCIParamAllConds    = [MCIParamAllConds,MCIParam];         
        end

        % remove Nan rows (because of 1) removing participants with short walking length; 2) not enough trials for parameter estimation)
        YoungParamAllConds      = removeNanRows(YoungParamAllConds);
        HealthyOldParamAllConds = removeNanRows(HealthyOldParamAllConds);
        [MCIParamAllConds MCInanIdx] = removeNanRows(MCIParamAllConds);

        YoungParamMean      = mean(YoungParamAllConds, 2);
        HealthyOldParamMean = mean(HealthyOldParamAllConds, 2);
        MCIParamMean        = mean(MCIParamAllConds, 2);
    
        f = figure('visible','off','Position', [100 100 500 500]);

        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)     

        colorForYoung = config.color_scheme_npg(3,:);        
        colorForHOld = config.color_scheme_npg(5,:);
        colorForMCI = config.color_scheme_npg(2,:);
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
        bp3 = boxplot(MCIParamMean, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 3);
        set(bp3,'linewidth',box_lineWidth);        

        %% Boxplot visual changes
        h = findobj(gca,'Tag','Box'); 
        %get the Young box
        patch(get(h(1),'XData'),get(h(1),'YData'),colorForMCI,'FaceAlpha',box_color_transparency);        
        %get the HelthyOld box
        patch(get(h(2),'XData'),get(h(2),'YData'),colorForHOld,'FaceAlpha',box_color_transparency);
        %get the MCI box
        patch(get(h(3),'XData'),get(h(3),'YData'),colorForYoung,'FaceAlpha',box_color_transparency);

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
                'MarkerFaceColor',colorForHOld, ...
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

        %% Scatter plot for data and mean (MCI negative)
        num_points = length(MCIParamMean(:,:));
        hold on
        x = 3*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, MCIParamMean(:,:), scatter_markerSize, ...
                'filled', ...
                'v', ...
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForMCINeg, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth);

        %add mean + errorbar
        mean_MCI = mean(MCIParamMean);
        sem_MCI = std(MCIParamMean)./sqrt(length(MCIParamMean));
        errorbar(3,mean_MCI,sem_MCI,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18); 
        hold on
        %add mean point
        scatter(3, mean_MCI, 3*scatter_markerSize, 'd',...
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
        alldata = [YoungParamMean;HealthyOldParamMean;MCIParamMean];
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
            'XTick'       , (1:3),... 
            'XLim'        , [0.5, 3.5],...
            'YLim'        , lowupYlim,...
            'YTick'       , yticks,...
            'XTickLabel'  , {'Young','Elderly','MCI-'},...
            'LineWidth'   , 1.0        );

        ylabel(ParamName(ParamIndx));

        % Extract pvalues from multiple comparison of group effect. Adding
        % it to the figure
        multicomp_result = multicomp_tab1{ParamIndx};
        % 1       2             3
        % MCI     HealthyOld    Young
        PvalueYoungvsHealthyOld = multicomp_result(3,6); % Young vs. HealthyOld see Two-way anova for details
        PvalueHealthyOldvsMCI = multicomp_result(1,6); % HealthyOld v.s. MCI vs.  see Two-way anova for details
        PvalueYoungvsMCI = multicomp_result(2,6); % Young vs. MCI see Two-way anova for details 

        title(strcat(['P12 = ',sprintf('%.2g',PvalueYoungvsHealthyOld)],...
              ['    P23 = ',sprintf('%.2g',PvalueHealthyOldvsMCI)],...
              ['    P13 = ',sprintf('%.2g',PvalueYoungvsMCI)]))
        
        %% Add significance bars
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
            H=adjustablesigstar(XXX,AllP_select);
        end

        %% export figure
        exportgraphics(f,config.ResultFolder+"/MergeCondsBox_"+ParamName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/MergeCondsBox_"+ParamName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end