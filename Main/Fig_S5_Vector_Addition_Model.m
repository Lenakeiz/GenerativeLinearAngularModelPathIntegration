%% Script to create output for Fig. S5 - comparisons between a vector addition model inspired by Harootonian et al., 2020
% Zilong Ji, UCL, 2022 zilong.ji@ucl.ac.uk
% Using a vector addition model inspired by Harootonian et al., 2020 we
% show that elderly participants forget the first leg of the oubound path,
% more than the second one. This is a feature that we implemented in our
% linear and angular generative model. The model proposed here is a vector
% addition version of the GLAMPI, featuring production error in linear and
% angular returns.

% Preparing the data
GLAMPI_PrepareBaseConfig;

% Preprocessing the data
GLAMPI_PreprocessData;

% Model fitting
config.ModelName        =   "beta1_beta2_g2_g3_sigma_nu";
config.ParamName        =   ["beta1", "beta2", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI;

% Preparing the output
config.ResultFolder     =   pwd + "/Output/FigS5/"+config.ModelName;
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Generating color scheme for our paper
ColorPattern; 

% Collecting information from output
AllYoungParams      =   YoungControls.Results.estimatedParams;
AllHealthyOldParams =   HealthyControls.Results.estimatedParams;
AllMCIPosParams     =   MCIPos.Results.estimatedParams;
AllMCINegParams     =   MCINeg.Results.estimatedParams;
AllMCIUnkParams     =   MCIUnk.Results.estimatedParams;

% BarScatter Plot for HealthyOld 
BoxPlotOfBeta1Beta2(AllHealthyOldParams, config, "Elderly");

% BarScatter Plot for MCI 
BoxPlotOfBeta1Beta2(AllMCIParams, config, "MCI");

% Final cleanup to leave workspace as the end of the Preprocessing stage.
% Remove if you want to take a look at the output data.
clearvars -except config YoungControls HealthyControls MCINeg MCIPos MCIUnk

%% Box Plot Of Fitted Parameters by avearaging three paramters from three conditions into one mean value
function BoxPlotOfBeta1Beta2(AllParams, config, name)
    
    numConds = 3; % environmental conditions
    ParamName = config.ParamName;

    ParamAllConds = [];
    for TRIAL_FILTER=1:numConds
        %% extract data for each condition
        Param     = AllParams{TRIAL_FILTER}(:,1);
        ParamAllConds = [ParamAllConds,Param];       
    end

    % remove Nan rows (because of 1) removing participants with short walking length; 2) not enough trials for parameter estimation)
    ParamAllConds = removeNanRows(ParamAllConds);

    Beta1 = mean(ParamAllConds, 2);

    ParamAllConds = []; 
    for TRIAL_FILTER=1:numConds
        %% extract data for each condition
        Param     = AllParams{TRIAL_FILTER}(:,2);
        ParamAllConds = [ParamAllConds,Param];       
    end

    % remove Nan rows (because of 1) removing participants with short walking length; 2) not enough trials for parameter estimation)
    ParamAllConds = removeNanRows(ParamAllConds);

    Beta2 = mean(ParamAllConds, 2);

    % Calculating difference between groups
    [h,pValue, stats] = ttest(Beta1, Beta2, "Tail", "left")

    f = figure('visible','off','Position', [100 100, 400, 600]);

    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12) 
    
    Color1 = [0.5,0.5,0.5];
    Color2 = [0.5,0.5,0.5];
        
    %% Scatter plot
    bar(1, mean(Beta1),'FaceColor', Color1, 'EdgeColor', 'k', 'BarWidth',0.4, 'FaceAlpha', 0.8);
    hold on
    bar(2, mean(Beta2),'FaceColor', Color2, 'EdgeColor', 'k', 'BarWidth',0.4, 'FaceAlpha', 0.8);
    
    %% Connecting points from same participant
    for i=1:length(Beta1)
        hold on
        plot([1,2], [Beta1(i), Beta2(i)], 'Color',[0,0,0,0.2], LineWidth=1);
    end
    
    %% Scatter plot
    hold on
    scatter(ones(length(Beta1),1), Beta1,100, 'MarkerEdgeColor','k', ...
        'MarkerFaceColor',Color1, 'MarkerEdgeAlpha',0.2, 'MarkerFaceAlpha',0.2)
    hold on
    scatter(2*ones(length(Beta2),1), Beta2,100, 'MarkerEdgeColor','k', ...
        'MarkerFaceColor',Color2, 'MarkerEdgeAlpha',0.2, 'MarkerFaceAlpha',0.2)
    
    %% Connect mean values
    hold on
    plot([1,2], [mean(Beta1), mean(Beta2)], 'Color',[0,0,0,1.0], LineWidth=2);
    
    hold on
    %% Reference line
    yl = yline(1);
    yl.Color = 'red';
    yl.LineWidth = 2;
    yl.LineStyle = '--';

    %% Figure post-processing
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.0 .0 .0], ...
        'YColor'      , [.0 .0 .0], ...
        'XTick'       , (1:2),... 
        'XTickLabel'  , {'Beta1','Beta2'},...
        'XLim'        , [0.5, 2.5],...
        'YLim'        , [0,1.8],... 
        'LineWidth'   , 1.0        );

    %% Add significance bars
    AllP = [pValue];
    Xval = [[1,2]];
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

    yline(1,Color='r',LineStyle='--',LineWidth=2);

    %% Export figure
    exportgraphics(f,config.ResultFolder+"/MergeCondsBox_beta1beta2.png",'Resolution',300);
    exportgraphics(f,config.ResultFolder+"/MergeCondsBox_beta1beta2.pdf",'Resolution',300, 'ContentType','vector');
end