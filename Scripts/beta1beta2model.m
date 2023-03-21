%% test beta1 beta2 model
%% Preparing the data
VAM_PrepareBaseConfig;

%% Preprocessing the data
VAM_PreprocessData;

%% Preparing the data and Slecting the Model

%choose a model from the model zoo
config.ModelName        =   "beta1_beta2_g2_g3_sigma_nu";
config.ParamName        =   ["beta1", "beta2", "g2", "g3", "sigma", "nu"];

% config.ModelName        =   "beta1_beta2_sigma_nu";
% config.ParamName        =   ["beta1", "beta2", "sigma", "nu"];

config.NumParams        =   length(config.ParamName);

% Run the model
VAM;

%% Model run completed, preparing the data for plotting figures
config.ResultFolder     =   pwd + "/Output/ModelFigures/"+config.ModelName+"/Young_HealthyOld_MCINeg";
%create storing folder for trajectory if not exist
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

%% Generating color scheme
ColorPattern; 

%% Getting Information from results:
AllYoungParams      =   YoungControls.Results.estimatedParams;
AllHealthyOldParams =   HealthyControls.Results.estimatedParams;
AllMCIPosParams     =   MCIPos.Results.estimatedParams;
AllMCINegParams     =   MCINeg.Results.estimatedParams;
AllMCIUnkParams     =   MCIUnk.Results.estimatedParams;

%%
[AllMCIParams, AllMCIParamsStatusIndex]= MergeMCI(AllMCIPosParams, AllMCINegParams, AllMCIUnkParams);

%% BarScatter Plot for Young 
BoxPlotOfBeta1Beta2(AllYoungParams, config, "Young");

%% BarScatter Plot for HealthyOld 
BoxPlotOfBeta1Beta2(AllHealthyOldParams, config, "Elderly");

%% BarScatter Plot for MCI 
BoxPlotOfBeta1Beta2(AllMCIParams, config, "MCI");

%% Box Plot Of Fitted Parameters by avearaging three paramters from three conditions into one mean value
function BoxPlotOfBeta1Beta2(AllParams, config, name)
    
    numConds = 3; %3 is the condition number
    ParamName = config.ParamName;

    ParamAllConds = []; %dimension are different, so separate from MCIParamAllConds
    for TRIAL_FILTER=1:numConds
        %% extract data
        Param     = AllParams{TRIAL_FILTER}(:,1);
        ParamAllConds = [ParamAllConds,Param];       
    end

    %remove Nan rows (nan coz of 1, removing participants with short walking length; 2, not enough trials for parameter estimation)
    ParamAllConds = removeNanRows(ParamAllConds);

    Beta1 = mean(ParamAllConds, 2);

    ParamAllConds = []; %dimension are different, so separate from MCIParamAllConds
    for TRIAL_FILTER=1:numConds
        %% extract data
        Param     = AllParams{TRIAL_FILTER}(:,2);
        ParamAllConds = [ParamAllConds,Param];       
    end

    %remove Nan rows (nan coz of 1, removing participants with short walking length; 2, not enough trials for parameter estimation)
    ParamAllConds = removeNanRows(ParamAllConds);

    Beta2 = mean(ParamAllConds, 2);

    signrank(Beta1,Beta2,'tail','left')

    f = figure('visible','on','Position', [100 100, 400, 600]);
    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12) 
    
    Color1 = [0.5,0.5,0.5];
    Color2 = [0.5,0.5,0.5];
    
    
    %bar scatter plot
    bar(1, mean(Beta1),'FaceColor', Color1, 'EdgeColor', 'k', 'BarWidth',0.4, 'FaceAlpha', 0.8);
    hold on
    bar(2, mean(Beta2),'FaceColor', Color2, 'EdgeColor', 'k', 'BarWidth',0.4, 'FaceAlpha', 0.8);
    
    %add horizontal connection lines
    for i=1:length(Beta1)
        hold on
        plot([1,2], [Beta1(i), Beta2(i)], 'Color',[0,0,0,0.2], LineWidth=1);
    end
    
    %scatter plot
    hold on
    scatter(ones(length(Beta1),1), Beta1,100, 'MarkerEdgeColor','k', ...
        'MarkerFaceColor',Color1, 'MarkerEdgeAlpha',0.2, 'MarkerFaceAlpha',0.2)
    hold on
    scatter(2*ones(length(Beta2),1), Beta2,100, 'MarkerEdgeColor','k', ...
        'MarkerFaceColor',Color2, 'MarkerEdgeAlpha',0.2, 'MarkerFaceAlpha',0.2)
    
    %link the mean value
    hold on
    plot([1,2], [mean(Beta1), mean(Beta2)], 'Color',[0,0,0,1.0], LineWidth=2);
    
    hold on
    %add a reference line y=1
    yl = yline(1);
    yl.Color = 'red';
    yl.LineWidth = 2;
    yl.LineStyle = '--';
    
    %[h,p, stats] = ttest(Beta1, Beta2, "Tail", "left")
    signrank(Beta1, Beta2, "Tail", "left")
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
    title(name)
    %% save figure
    exportgraphics(f,config.ResultFolder+"/ZMergeCondsBox_beta1beta2.png",'Resolution',300);
    exportgraphics(f,config.ResultFolder+"/ZMergeCondsBox_beta1beta2.pdf",'Resolution',300, 'ContentType','vector');
end