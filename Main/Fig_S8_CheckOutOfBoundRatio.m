%% Script to check the out-of-bound ratio in different groups of data
% Zilong Ji, UCL, 2022 zilong.ji@ucl.ac.uk

% Preparing the data
GLAMPI_PrepareBaseConfig

% Preprocessing the data
GLAMPI_PreprocessData

%%

%get the out-of-bound ratio for each participant under each condition
OoBratioYoung = getOoBRatio(YoungControls);
OoBratioHO = getOoBRatio(HealthyControls);
OoBratioMCIUnk = getOoBRatio(MCIUnk);
OoBratioMCIPos = getOoBRatio(MCIPos);
OoBratioMCINeg = getOoBRatio(MCINeg);

config.ResultFolder = pwd + "/Output/FigS8";
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

%% perform a two-way anova
[anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnova_OoB(OoBratioYoung, OoBratioHO, OoBratioMCIUnk, OoBratioMCIPos, OoBratioMCINeg, config);

%% plot
BoxPlotOfMergeCondition(OoBratioYoung, OoBratioHO, OoBratioMCIUnk, OoBratioMCIPos, OoBratioMCINeg, multicomp_tab1, config)
%%
function OoBRatio_all = getOoBRatio(Group)
    num = length(Group.CondTable);
    
    bpids = Group.BadPptIdxs;
    OoBRatio_all = zeros(num-length(bpids),3);

    ind = 1;
    num_oob = 0;
    num_total = 0;
    for i = 1:num
        if ismember(i, bpids)
            continue;
        end
        participant_i = Group.CondTable{i};
        BadExecution_i = Group.Reconstructed{i}.BadExecution;

        conditions_i = participant_i.Condition;
        oob_i = participant_i.OutOfBound;
        for cond = 1:3
            oob_i_c = oob_i(conditions_i == cond);
            bexc_i_c = BadExecution_i(conditions_i == cond);

            %select the first 9 trials for young
            if length(oob_i_c)>10
                oob_i_c = oob_i_c(1:9);
                bexc_i_c = bexc_i_c(1:9);
            end
            %filter oob_i_c with BadExecution=0
            oob_i_c = oob_i_c(bexc_i_c==0);
            %each participant, under condition cond, calculate the OoB
            %ratio under that condition
            oob_i_c_ratio = sum(oob_i_c==1)/length(oob_i_c);

            OoBRatio_all(ind, cond) = oob_i_c_ratio;

            num_oob = num_oob + sum(oob_i_c==1);
            num_total = num_total + length(oob_i_c);
        end
        ind = ind +1;
    end
    num_oob
    num_total
    num_oob/num_total
end

%%
function BoxPlotOfMergeCondition(Young, HO, MCIUnk, MCIPos, MCINeg, multicomp_tab1, config)
    
    YoungMean = mean(Young,2);
    HOMean = mean(HO,2);
    MCIUnkMean = mean(MCIUnk,2);
    MCIPosMean = mean(MCIPos,2);
    MCINegMean = mean(MCINeg,2);

    f = figure('visible','off','Position', [100 100 600 300]);

    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)     

    colorForYoung = config.color_scheme_npg(3,:);        
    colorForHOld = config.color_scheme_npg(5,:);
    colorForMCIUnk = config.color_scheme_npg(2,:);
    colorForMCINeg = config.color_scheme_npg(4,:);
    colorForMCIPos = config.color_scheme_npg(6,:);

    % parameters set for controlling visual output
    whisker_value               =   1.5;
    box_lineWidth               =   0.3;
    box_widths_value            =   0.4;
    box_color_transparency      =   0.5;
    median_lineWidth            =   2;
    median_color                =   'k';
    scatter_jitter_value        =   0.2;
    scatter_markerSize          =   40;
    scatter_marker_edgeColor    =   'k';
    scatter_marker_edgeWidth    =   0.5;
    scatter_color_transparency  =   0.7;

    hold on
    %% Boxplot for each column in Young
    bp1 = boxplot(YoungMean, ...
                'Whisker',whisker_value, ...
                'symbol','', ... %symbol ='' making outlier invisible
                'Color','k', ...
                'Notch','on', ...
                'widths',box_widths_value,...
                'positions', 1);
    set(bp1,'linewidth',box_lineWidth);

    hold on
    %% Boxplot for each column in Healthy Elderly
    bp2 = boxplot(HOMean, ...
                'Whisker',whisker_value, ...
                'symbol','', ... %symbol ='' making outlier invisible
                'Color','k', ...
                'Notch','on', ...
                'widths',box_widths_value,...
                'positions', 2);
    set(bp2,'linewidth',box_lineWidth);

    hold on
    %% Boxplot for each column in MCI Unk
    bp3 = boxplot(MCIUnkMean, ...
                'Whisker',whisker_value, ...
                'symbol','', ... %symbol ='' making outlier invisible
                'Color','k', ...
                'Notch','on', ...
                'widths',box_widths_value,...
                'positions', 3);
    set(bp3,'linewidth',box_lineWidth);        

    hold on
    %% Boxplot for each column in MCI Pos
    bp4 = boxplot(MCIPosMean, ...
                'Whisker',whisker_value, ...
                'symbol','', ... %symbol ='' making outlier invisible
                'Color','k', ...
                'Notch','on', ...
                'widths',box_widths_value,...
                'positions', 4);
    set(bp4,'linewidth',box_lineWidth);        

    hold on
    %% Boxplot for each column in MCI Neg
    bp5 = boxplot(MCINegMean, ...
                'Whisker',whisker_value, ...
                'symbol','', ... %symbol ='' making outlier invisible
                'Color','k', ...
                'Notch','on', ...
                'widths',box_widths_value,...
                'positions', 5);
    set(bp5,'linewidth',box_lineWidth); 

    %% Boxplot visual changes
    h = findobj(gca,'Tag','Box'); 
    patch(get(h(1),'XData'),get(h(1),'YData'),colorForMCIPos,'FaceAlpha',box_color_transparency);        
    patch(get(h(2),'XData'),get(h(2),'YData'),colorForMCINeg,'FaceAlpha',box_color_transparency);
    patch(get(h(3),'XData'),get(h(3),'YData'),colorForMCIUnk,'FaceAlpha',box_color_transparency);
    patch(get(h(4),'XData'),get(h(4),'YData'),colorForHOld,'FaceAlpha',box_color_transparency);
    patch(get(h(5),'XData'),get(h(5),'YData'),colorForYoung,'FaceAlpha',box_color_transparency);
    
    %% Median visual change
    h=findobj(gca,'tag','Median');
    for i = 1:length(h)
        h(i).LineWidth = median_lineWidth;
        h(i).Color = median_color;
    end

    %% Scatter plot for data and mean (Young)
    num_points = size(YoungMean,1);
    hold on
    x = 1*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
    scatter(x, YoungMean, scatter_markerSize, ...
            'filled', ...
            'o', ... %marker shape
            'MarkerEdgeColor',scatter_marker_edgeColor, ...
            'MarkerFaceColor',colorForYoung, ...
            'MarkerFaceAlpha',scatter_color_transparency,...
            'LineWidth',scatter_marker_edgeWidth);   

    %% Scatter plot for data and mean (Healthy Elderly)
    num_points = length(HOMean);
    hold on
    x = 2*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
    scatter(x, HOMean, scatter_markerSize, ...
            'filled', ...
            'o', ... %marker shape
            'MarkerEdgeColor',scatter_marker_edgeColor, ...
            'MarkerFaceColor',colorForHOld, ...
            'MarkerFaceAlpha',scatter_color_transparency,...
            'LineWidth',scatter_marker_edgeWidth);   

    %% Scatter plot for data and mean (MCI unknown)
    num_points = length(MCIUnkMean);
    hold on
    x = 3*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); 
    scatter(x, MCIUnkMean, scatter_markerSize, ...
            'filled', ...
            'o', ...
            'MarkerEdgeColor',scatter_marker_edgeColor, ...
            'MarkerFaceColor',colorForMCIUnk, ...
            'MarkerFaceAlpha',scatter_color_transparency,...
            'LineWidth',scatter_marker_edgeWidth);


    %% Scatter plot for data and mean (MCI unknown)
    num_points = length(MCIPosMean);
    hold on
    x = 4*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); 
    scatter(x, MCIPosMean, scatter_markerSize, ...
            'filled', ...
            'o', ...
            'MarkerEdgeColor',scatter_marker_edgeColor, ...
            'MarkerFaceColor',colorForMCIPos, ...
            'MarkerFaceAlpha',scatter_color_transparency,...
            'LineWidth',scatter_marker_edgeWidth);

    %% Scatter plot for data and mean (MCI unknown)
    num_points = length(MCINegMean);
    hold on
    x = 5*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); 
    scatter(x, MCINegMean, scatter_markerSize, ...
            'filled', ...
            'o', ...
            'MarkerEdgeColor',scatter_marker_edgeColor, ...
            'MarkerFaceColor',colorForMCINeg, ...
            'MarkerFaceAlpha',scatter_color_transparency,...
            'LineWidth',scatter_marker_edgeWidth);

    %% Figure post-processing

    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.0 .0 .0], ...
        'YColor'      , [.0 .0 .0], ...
        'XTick'       , (1:5),... 
        'XLim'        , [0.5, 5.5],...
        'XTickLabel'  , {'Young','Elderly','MCIUnk', 'MCIPos', 'MCINeg'},...
        'FontSize'    , 12,...
        'LineWidth'   , 1.0        );

    ylabel('Out-of-boundary Ratio');
    xtickangle(45);

    % Extract pvalues from multiple comparison of group effect. Adding
    % it to the figure
    multicomp_result = multicomp_tab1;
    % 1       2       3        4        5
    % Young   HO    MCIUnk   MCIPos  MCINeg
    PvalueYoungvsHealthyOld = multicomp_result(1,6); % Young vs. HealthyOld see Two-way anova for details
    PvalueYoungvsMCIUnk= multicomp_result(2,6); % HealthyOld v.s. MCI vs.  see Two-way anova for details
    PvalueYoungvsMCIPos = multicomp_result(3,6); % Young vs. MCI see Two-way anova for details 
    PvalueYoungvsMCINeg = multicomp_result(4,6); % Young vs. MCI see Two-way anova for details 

        
    %% Add significance bars 
    AllP = [PvalueYoungvsHealthyOld,PvalueYoungvsMCIUnk,PvalueYoungvsMCIPos, PvalueYoungvsMCINeg];
    Xval = [[1,2];[1,3];[1,4];[1,5]];
    %AllP = PvalueYoungvsHealthyOld;
    %Xval = [1,2];
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
        exportgraphics(f,config.ResultFolder+"/MergeCondsBox.png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/MergeCondsBox.pdf",'Resolution',300, 'ContentType','vector');

end