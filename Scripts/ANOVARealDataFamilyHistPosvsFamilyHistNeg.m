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
config.NumParams                                        = 1000;
config.ModelName                                        = "beta_g2_g3_sigma_nu";

resultfolder = pwd + "/Output/DataFigures/FamilyHistPos_FamilyHistNeg";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% Pos data
FamilyHistPos   = TransformPaths(FamilyHistPos);%transform data
FamilyHistPos   = CalculateTrackingPath(FamilyHistPos, config);
ManuallyScoringFamilyHistPos;
%%
[~, FHPosX, FHPosDX, FHPosTheta, FHPosDistErr, FHPosAngleErr, FHPosFlagOoB, ~] = getResultsAllConditions(FamilyHistPos, config);

%% Model fitting for Neg data
FamilyHistNeg   = TransformPaths(FamilyHistNeg);%transform data
FamilyHistNeg   = CalculateTrackingPath(FamilyHistNeg, config);
ManuallyScoringFamilyHistNeg;
%%
[~, FHNegX, FHNegDX, FHNegTheta, FHNegDistErr, FHNegAngleErr, FHNegFlagOoB, ~] = getResultsAllConditions(FamilyHistNeg, config);

%% Setting colors for using in plots
ColorPattern; 

%% Error plot
GenderPos = FamilyHistPos.Gender;
GenderNeg = FamilyHistNeg.Gender;
ErrPlot(FHPosDistErr, FHNegDistErr, GenderPos, GenderNeg, 'dist', config)
ErrPlot(FHPosAngleErr, FHNegAngleErr, GenderPos, GenderNeg, 'angle', config)

%% Scatter Error Plot in Male 
ScatterErrPlotinMale(FHPosAngleErr, FHNegAngleErr, GenderPos, GenderNeg, FHPosFlagOoB, FHNegFlagOoB, 'angle', config)

%%
function ScatterErrPlotinMale(FHPosErr, FHNegErr, GenderPos, GenderNeg, FHPosFlagOoB, FHNegFlagOoB, type, config)
    %% Pos 
    PosmaleIdx = GenderPos==1;

    PosMale = zeros(sum(PosmaleIdx),3);
    ScatterPosMale_OoB = cell(1,3);
    ScatterPosMale_InB = cell(1,3);

    for cond=1:3
        ScatterCond_OoB = [];
        ScatterCond_InB = [];
        PosMale_cond_i = FHPosErr{cond}(PosmaleIdx);
        PosMale_cond_i_OoB = FHPosFlagOoB{cond}(PosmaleIdx);
        for id=1:length(PosMale_cond_i)
            err = PosMale_cond_i{id};
            err = cell2mat(err);
            flagOoB = logical(PosMale_cond_i_OoB{id});

            InB_err = err(~flagOoB);
            meanerr = mean(InB_err, 'omitnan');
            PosMale(id,cond) = meanerr;

            OoB_err = err(flagOoB);
            ScatterCond_OoB = [ScatterCond_OoB, OoB_err];
            ScatterCond_InB = [ScatterCond_InB, InB_err];
        end 

        ScatterPosMale_OoB{cond} = ScatterCond_OoB;
        ScatterPosMale_InB{cond} = ScatterCond_InB;
    end

    %swap last two column
    v = PosMale(:,2); PosMale(:,2)=PosMale(:,3); PosMale(:,3)=v;


    %% Neg
    NegmaleIdx = GenderNeg==1;

    NegMale = zeros(sum(NegmaleIdx),3);
    ScatterNegMale_OoB = cell(1,3);
    ScatterNegMale_InB = cell(1,3);

    for cond=1:3
        ScatterCond_OoB = [];
        ScatterCond_InB = [];
        NegMale_cond_i = FHNegErr{cond}(NegmaleIdx);
        NegMale_cond_i_OoB = FHNegFlagOoB{cond}(NegmaleIdx);
        for id=1:length(NegMale_cond_i)
            err = NegMale_cond_i{id};
            err = cell2mat(err);
            flagOoB = logical(NegMale_cond_i_OoB{id});

            InB_err = err(~flagOoB);
            meanerr = mean(InB_err, 'omitnan');
            NegMale(id,cond) = meanerr;

            OoB_err = err(flagOoB);
            ScatterCond_OoB = [ScatterCond_OoB, OoB_err];
            ScatterCond_InB = [ScatterCond_InB, InB_err];
        end 

        ScatterNegMale_OoB{cond} = ScatterCond_OoB;
        ScatterNegMale_InB{cond} = ScatterCond_InB;
    end

    %swap last two column
    v = NegMale(:,2); NegMale(:,2)=NegMale(:,3); NegMale(:,3)=v;

    %% set figure info
    f = figure('visible','on','Position', [100 100 600 600]);
    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)   

    colorForPos = config.color_scheme_npg(5,:);
    colorForNeg = config.color_scheme_npg(2,:);

    %% male
    meanPosMale = mean(PosMale, 1, 'omitnan');
    semPosMale = std(PosMale, 1, 'omitnan')/sqrt(size(PosMale,1));
    meanNegMale = mean(NegMale, 1, 'omitnan');
    semNegMale = std(NegMale, 1, 'omitnan')/sqrt(size(NegMale,1));
    errorbar([1,2,3], meanPosMale, semPosMale, '-o', 'MarkerSize',10, 'Color', colorForPos, 'LineWidth',3);
    hold on
    errorbar([1,2,3], meanNegMale, semNegMale, '-o', 'MarkerSize',10, 'Color', colorForNeg, 'LineWidth',3);

    if type=="dist"
        ylabel('Distance error (meters)')
        ylim = [0,1.2];
        ytick = [0,0.3,0.6,0.9];
    elseif type=="angle"
        ylabel('Angular error (rads)')
        ylim = [-30,30];
        ytick=[-30,0,30];
    else
        %
    end

    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , (1:3),... 
        'YTick'       , ytick,... 
        'XTickLabel'  , {'No Change','No Optical Flow', 'No Distal Cue'},...
        'XLim'        , [0.5, 3.5],...
        'YLim'        , ylim,...
        'LineWidth'   , .5        );

    title('Male')
    legend({'FHPos','FHNeg'})
end

%%
function ErrPlot(FHPosErr, FHNegErr, GenderPos, GenderNeg, type, config)
    %% Pos 
    PosmaleIdx = GenderPos==1;
    PosfemaleIdx = GenderPos==2;

    PosMale = zeros(sum(PosmaleIdx),3);
    PosFemale = zeros(sum(PosfemaleIdx),3);

    for cond=1:3
        PosMale_cond_i = FHPosErr{cond}(PosmaleIdx);
        for id=1:length(PosMale_cond_i)
            err = PosMale_cond_i{id};
            err = cell2mat(err);
            meanerr = mean(err, 'omitnan');
            PosMale(id,cond) = meanerr;
        end

        PosFemale_cond_i = FHPosErr{cond}(PosfemaleIdx);
        for id=1:length(PosFemale_cond_i)
            err = PosFemale_cond_i{id};
            err = cell2mat(err);
            meanerr = mean(err, 'omitnan');
            PosFemale(id,cond) = meanerr;
        end   
    end

    %swap last two column
    v = PosMale(:,2); PosMale(:,2)=PosMale(:,3); PosMale(:,3)=v;
    %swap last two column
    v = PosFemale(:,2); PosFemale(:,2)=PosFemale(:,3); PosFemale(:,3)=v;

    
    %% Neg
    NegmaleIdx = GenderNeg==1;
    NegfemaleIdx = GenderNeg==2;

    NegMale = zeros(sum(NegmaleIdx),3);
    NegFemale = zeros(sum(NegfemaleIdx),3);

    for cond=1:3
        NegMale_cond_i = FHNegErr{cond}(NegmaleIdx);
        for id=1:length(NegMale_cond_i)
            err = NegMale_cond_i{id};
            err = cell2mat(err);
            meanerr = mean(err, 'omitnan');
            NegMale(id,cond) = meanerr;
        end

        NegFemale_cond_i = FHNegErr{cond}(NegfemaleIdx);
        for id=1:length(NegFemale_cond_i)
            err = NegFemale_cond_i{id};
            err = cell2mat(err);
            meanerr = mean(err, 'omitnan');
            NegFemale(id,cond) = meanerr;
        end 
    end

    %swap last two column
    v = NegMale(:,2); NegMale(:,2)=NegMale(:,3); NegMale(:,3)=v;
    %swap last two column
    v = NegFemale(:,2); NegFemale(:,2)=NegFemale(:,3); NegFemale(:,3)=v;


    %% set figure info
    f = figure('visible','on','Position', [100 100 600 300]);
    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)   

    colorForPos = config.color_scheme_npg(5,:);
    colorForNeg = config.color_scheme_npg(2,:);

    %% male
    subplot(1,2,1);
    meanPosMale = mean(PosMale, 1, 'omitnan');
    semPosMale = std(PosMale, 1, 'omitnan')/sqrt(size(PosMale,1));
    meanNegMale = mean(NegMale, 1, 'omitnan');
    semNegMale = std(NegMale, 1, 'omitnan')/sqrt(size(NegMale,1));
    errorbar([1,2,3], meanPosMale, semPosMale, '-o', 'MarkerSize',10, 'Color', colorForPos, 'LineWidth',3);
    hold on
    errorbar([1,2,3], meanNegMale, semNegMale, '-o', 'MarkerSize',10, 'Color', colorForNeg, 'LineWidth',3);

    if type=="dist"
        ylabel('Distance error (meters)')
        ylim = [0,1.2];
        ytick = [0,0.3,0.6,0.9];
    elseif type=="angle"
        ylabel('Angular error (rads)')
        %ylim = [-30,30];
        %ytick=[0,10,20];
        ylim = [-10,20];
        ytick=[-10,0,10,20];
    else
        %
    end

    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , (1:3),... 
        'YTick'       , ytick,... 
        'XTickLabel'  , {'No Change','No Optical Flow', 'No Distal Cue'},...
        'XLim'        , [0.5, 3.5],...
        'YLim'        , ylim,...
        'LineWidth'   , .5        );

    title('Male')
    legend({'FHPos','FHNeg'})

    %% female
    subplot(1,2,2);
    meanPosFemale = mean(PosFemale, 1, 'omitnan');
    semPosFemale = std(PosFemale, 1, 'omitnan')/sqrt(size(PosFemale,1));
    meanNegFemale = mean(NegFemale, 1, 'omitnan');
    semNegFemale = std(NegFemale, 1, 'omitnan')/sqrt(size(NegFemale,1));
    errorbar([1,2,3], meanPosFemale, semPosFemale, '-o', 'MarkerSize',10, 'Color', colorForPos, 'LineWidth',3);
    hold on
    errorbar([1,2,3], meanNegFemale, semNegFemale, '-o', 'MarkerSize',10, 'Color', colorForNeg, 'LineWidth',3);
  
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , (1:3),... 
        'YTick'       , ytick,... 
        'XTickLabel'  , {'No Change', 'No Optical Flow', 'No Distal Cue'},...
        'XLim'        , [0.5, 3.5],...
        'YLim'        , ylim,...        
        'LineWidth'   , .5        );

    title('Female')
    legend({'FHPos','FHNeg'})

end

%% TwowayAnova and BarScatter Plot
% [anova_tab,multicomp_tab1,~, ~] = TwowayAnova_CocoData(AllFamilyHistPosParams, AllFamilyHistNegParams, config);
% BoxPlotOfFittedParam(AllFamilyHistPosParams, AllFamilyHistNegParams, anova_tab, config);
% BoxPlotOfFittedParamMergeCondition(AllFamilyHistPosParams, AllFamilyHistNegParams, multicomp_tab1, config)

%% ThreewayAnova
% config.Gender_Pos = FamilyHistPos.Gender;
% config.Gender_Neg = FamilyHistNeg.Gender;
% [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab3, multicomp_tab12] = ThreewayAnova_CocoData(AllFamilyHistPosParams, AllFamilyHistNegParams, config);

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