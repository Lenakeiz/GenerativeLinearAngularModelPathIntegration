%% Script to create output for Fig. 3 - visualization of the GLAMPI model applied to different groups
% Zilong Ji, UCL, 2022 zilong.ji@ucl.ac.uk
% Andrea Castegnaro, UCL, 2022 uceeaca@ucl.ac.uk
% Visualization include, visualizing mental points over real data. 
% Comparison between effect ratio between mental/real leg 1 and leg 2.
% Comparison between real and mental turn between leg 1 and leg 2.
% Regression to the distance mean effect as indicated by m3.
% Regression to the angular mean effect indicated by g3.
% Plots the distribution of each single parameter.

% Preparing the data
GLAMPI_PrepareBaseConfig;

% Preprocessing the data
GLAMPI_PreprocessData;

% Model fitting
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI;

% Preparing the output
config.ResultFolder     =   pwd + "/Output/Fig3/"+config.ModelName+"/VisualizingMentalPhysicalTrajectory";
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Generating color scheme for our paper
ColorPattern; 

% Collecting information from output
AllYoungResults         =   YoungControls.Results;
AllHealthyOldResults    =   HealthyControls.Results;
AllMCIPosResults        =   MCIPos.Results;
AllMCINegResults        =   MCINeg.Results;
AllMCIUnkResults        =   MCIUnk.Results;

%% Starting visualization
% Can select a different groups by selecting different names (creates Fig S6, S7)
% name = "Young";
name = "HealthyOld";
% name = "MCIMerged";

%% Visualizing trials (Healthy elderly only)
id = 13; % participant Id
trials = [1,2,3,4,5]; % trial numbers

for tN=1:length(trials)
    trial_id = trials(tN);
    [phy_p3,men_p3] = VisualizeMenPhyTraj(AllHealthyOldResults, id, 1, trial_id, true, config); %true means plot the figure
end

% Visualize some trials
id = 28; % participant Id
trials = [1,2,3,4,5]; % trial numbers

for tN=1:length(trials)
    trial_id = trials(tN);
    [phy_p3,men_p3] = VisualizeMenPhyTraj(AllHealthyOldResults, id, 1, trial_id, true, config); %true means plot the figure
end

%% Visualizing group scatter plot of real returned point and mental return point

if name == "Young"
    [Physical_Pos, Mental_Pos, L1Diff, L2Diff, L1Ratio, L2Ratio, T2, T2_prime, Alpha, T3_prime] = extractAllFinalPoints(AllYoungResults, config);
elseif name=="HealthyOld"
    [Physical_Pos, Mental_Pos, L1Diff, L2Diff,  L1Ratio, L2Ratio, T2, T2_prime, Alpha, T3_prime] = extractAllFinalPoints(AllHealthyOldResults, config);
elseif name=="MCIMerged"
    [Physical_Pos1, Mental_Pos1, L1Diff1, L2Diff1, L1Ratio1, L2Ratio1, T2_1, T2_prime1, Alpha1, T3_prime1] = extractAllFinalPoints(AllMCIPosResults, config);
    [Physical_Pos2, Mental_Pos2, L1Diff2, L2Diff2, L1Ratio2, L2Ratio2, T2_2, T2_prime2, Alpha2, T3_prime2] = extractAllFinalPoints(AllMCINegResults, config);
    [Physical_Pos3, Mental_Pos3, L1Diff3, L2Diff3, L1Ratio3, L2Ratio3, T2_3, T2_prime3, Alpha3, T3_prime3] = extractAllFinalPoints(AllMCIUnkResults, config);
    Physical_Pos = [Physical_Pos1;Physical_Pos2;Physical_Pos3];
    Mental_Pos = [Mental_Pos1;Mental_Pos2;Mental_Pos3];
else
    error("Choose correct name!")
end

% Plotting physical and mental points
f = figure('visible','off','Position', [100 100 500 500]);

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12) 

MenColor = config.color_scheme_npg(1,:);
PhyColor = config.color_scheme_npg(2,:);

% Plotting real positions
scatter(Physical_Pos(:,1), Physical_Pos(:,2), 'MarkerEdgeColor',PhyColor, ...
    'MarkerFaceColor',PhyColor,'MarkerEdgeAlpha',0.2, 'MarkerFaceAlpha',0.4);
hold on
% Plotting mental positions
scatter(Mental_Pos(:,1), Mental_Pos(:,2), 'MarkerEdgeColor',MenColor, ...
    'MarkerFaceColor', MenColor, 'MarkerEdgeAlpha',0.2, 'MarkerFaceAlpha',0.4);

legend('Actual return points', 'Generated return points');

% Figure post-processing
set(gca, ...
'Box'         , 'off'       , ...
'TickDir'     , 'out'       , ...
'TickLength'  , [.01 .01]   , ...
'XColor'      , [.1 .1 .1]  , ...
'YColor'      , [.1 .1 .1]  , ...
'XLim'        , [-3.5, 6.5] ,...
'YLim'        , [-3.5, 6.5] ,... 
'XTick'       , [-3,0,3,6]    ,...
'YTick'       , [-3,0,3,6]    ,...
'LineWidth'   , .5          ,...
'XAxisLocation', 'origin'   ,...
'YAxisLocation', 'origin');

% Export figure
exportgraphics(f,config.ResultFolder+"/"+name+".png",'Resolution',300);
exportgraphics(f,config.ResultFolder+"/"+name+".pdf",'Resolution',300, 'ContentType','vector');

%% Plotting regression to the mean distance effect (m3)
f = figure('visible','off','Position', [100 100 500 500]);

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12) 

MenColor = config.color_scheme_npg(1,:);
PhyColor = config.color_scheme_npg(2,:);

Physical_Dist = sqrt(sum(Physical_Pos.^2,2));
Mental_Dist = sqrt(sum(Mental_Pos.^2,2));

scatter(Physical_Dist, Mental_Dist,50, 'MarkerEdgeColor','k', ...
    'MarkerFaceColor','k', 'MarkerEdgeAlpha',0.1, 'MarkerFaceAlpha',0.1)
sh = scatterhist(Physical_Dist,Mental_Dist, 'Kernel','overlay','Location','NorthEast', 'Color','k', 'Direction','out')
x = get(gca,'children');
set(x,'markerfacecolor',[0.9,0.9,0.9])

hold on
c = Physical_Dist\Mental_Dist;
xc = linspace(0,6,10);
plot(xc,c*xc,'r--', LineWidth=2);

hold on

x = linspace(0,6,10); y = x; 
plot(x,y,'k--', LineWidth=2);

xlabel("Actual distance to origin (meters)")
ylabel("Mental distance to origin (meters)");

set(gca, ...
'Box'         , 'off'       , ...
'TickDir'     , 'out'       , ...
'TickLength'  , [.01 .01]   , ...
'XColor'      , [.1 .1 .1]  , ...
'YColor'      , [.1 .1 .1]  , ...
'XLim'        , [0, 6] ,...
'YLim'        , [0, 6] ,... 
'XTick'       , [0,2,4,6]    ,...
'YTick'       , [0,2,4,6]    ,...
'LineWidth'   , .5          ,...
'XAxisLocation', 'origin'   ,...
'YAxisLocation', 'origin');

exportgraphics(f,config.ResultFolder+"/"+name+"_distribution.png",'Resolution',300);
exportgraphics(f,config.ResultFolder+"/"+name+"_distribution.pdf",'Resolution',300, 'ContentType','vector');

%% Bar plot of comparing forgetting amount on leg 1 versus leg 2 (beta)
noNan_L1Ratio = L1Ratio(~isnan(L1Ratio));
noNan_L2Ratio = L2Ratio(~isnan(L2Ratio));

f = figure('visible','off','Position', [100 100, 500, 400]);
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12) 

Color1 = [0.5,0.5,0.5];
Color2 = [0.5,0.5,0.5];

% Bar plots
bar(1, mean(noNan_L1Ratio),'FaceColor', Color1, 'EdgeColor', 'k', 'BarWidth',0.4, 'FaceAlpha', 0.8);
hold on
bar(2, mean(noNan_L2Ratio),'FaceColor', Color2, 'EdgeColor', 'k', 'BarWidth',0.4, 'FaceAlpha', 0.8);

% Adding connection lines
for i=1:length(noNan_L1Ratio)
    hold on
    plot([1,2], [noNan_L1Ratio(i), noNan_L2Ratio(i)], 'Color',[0,0,0,0.2], LineWidth=1);
end

% Data scatter points
hold on
scatter(ones(length(noNan_L1Ratio),1), noNan_L1Ratio,100, 'MarkerEdgeColor','k', ...
    'MarkerFaceColor',Color1, 'MarkerEdgeAlpha',0.2, 'MarkerFaceAlpha',0.2)
hold on
scatter(2*ones(length(noNan_L2Ratio),1), noNan_L2Ratio,100, 'MarkerEdgeColor','k', ...
    'MarkerFaceColor',Color2, 'MarkerEdgeAlpha',0.2, 'MarkerFaceAlpha',0.2)

% Link mean values
hold on
plot([1,2], [mean(noNan_L1Ratio), mean(noNan_L2Ratio)], 'Color',[0,0,0,1.0], LineWidth=2);

hold on

% Reference line
yl = yline(1);
yl.Color = 'red';
yl.LineWidth = 2;
yl.LineStyle = '--';

[h,p, ci, stats] = ttest(noNan_L1Ratio, noNan_L2Ratio, "Tail", "left")

set(gca, ...
'Box'         , 'off'       , ...
'TickDir'     , 'out'       , ...
'TickLength'  , [.01 .01]   , ...
'XColor'      , [.1 .1 .1]  , ...
'YColor'      , [.1 .1 .1]  , ...
'XLim'        , [0.5, 2.5] ,...
'YLim'        , [0,1.6],...
'XTick'       , [1,2]    ,...
'YTick'       , [0,0.5,1.0,1.5] ,...
'XTickLabel'   , {'Leg 1', 'Leg 2'},...
'LineWidth'   , .5          ,...
'XAxisLocation', 'origin'   ,...
'YAxisLocation', 'origin');

ylabel("Encoded distance / actual distance");

exportgraphics(f,config.ResultFolder+"/"+name+"_legratio.png",'Resolution',300);
exportgraphics(f,config.ResultFolder+"/"+name+"_legratio.pdf",'Resolution',300, 'ContentType','vector');

%% Plotting real turn vs model generated turn between first and second segment (g2)
noNan_T2 = T2(~isnan(T2));
noNan_T2_prime = T2_prime(~isnan(T2_prime));

f = figure('visible','off','Position', [100 100, 500, 400]);
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12) 

MenColor = config.color_scheme_npg(1,:);
PhyColor = config.color_scheme_npg(2,:);

% Bar plots
bar(1, mean(noNan_T2),'FaceColor', PhyColor, 'EdgeColor', 'k', 'BarWidth',0.4, 'FaceAlpha', 0.8);
hold on
bar(2, mean(noNan_T2_prime),'FaceColor', MenColor, 'EdgeColor', 'k', 'BarWidth',0.4, 'FaceAlpha', 0.8);

hold on
% Mean data plotting
mean_noNan_T2 = mean(noNan_T2);
sem_noNan_T2 = std(noNan_T2)./sqrt(length(noNan_T2));
errorbar(1,mean_noNan_T2,sem_noNan_T2,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18); 
hold on
scatter(1, mean(noNan_T2), 200, 'MarkerEdgeColor','k', 'MarkerFaceColor',PhyColor, 'MarkerFaceAlpha',0.5)

hold on
mean_noNan_T2_prime = mean(noNan_T2_prime);
sem_noNan_T2_prime = std(noNan_T2_prime)./sqrt(length(noNan_T2_prime));
errorbar(2,mean_noNan_T2_prime,sem_noNan_T2_prime,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18); 
hold on
scatter(2, mean(noNan_T2_prime), 200, 'MarkerEdgeColor','k', 'MarkerFaceColor',MenColor, 'MarkerFaceAlpha',0.5)
  
hold on
plot([1,2], [mean(noNan_T2), mean(noNan_T2_prime)], 'Color',[0,0,0,1.0], LineWidth=2);

set(gca, ...
'Box'         , 'off'       , ...
'TickDir'     , 'out'       , ...
'TickLength'  , [.01 .01]   , ...
'XColor'      , [.1 .1 .1]  , ...
'YColor'      , [.1 .1 .1]  , ...
'XLim'        , [0.5, 2.5] ,...
'XTick'       , [1,2]    ,...
'XTickLabel'   , {'Actual turning \theta_2', 'Generated turning \theta_2^\prime'},...
'LineWidth'   , .5          ,...
'XAxisLocation', 'origin'   ,...
'YAxisLocation', 'origin');

ylabel('Turning angle (radians)')

[h,p, stats] = ttest(noNan_T2, noNan_T2_prime, "Tail", "right")

exportgraphics(f,config.ResultFolder+"/"+name+"_turn2.png",'Resolution',300);
exportgraphics(f,config.ResultFolder+"/"+name+"_turn2.pdf",'Resolution',300, 'ContentType','vector');

%% regression to the angular mean effect (g3)
f = figure('visible','off','Position', [100 100, 500, 400]);
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12) 

cmap = config.color_scheme_npg([3,4,6,7,9,10],:);

CMs = [];
for i=1:length(Alpha)
    if ~isempty(Alpha{i})
        mapidx = mod(i, length(cmap))+1;
        scatter(Alpha{i}, T3_prime{i}, 'MarkerEdgeColor',cmap(mapidx,:), 'MarkerFaceColor', cmap(mapidx,:), ...
            'MarkerFaceAlpha', 0., 'MarkerEdgeAlpha', 0.);
        CMs = [CMs;cmap(mapidx,:)];
    end
    hold on
end

ls = lsline;
X = [];
Y = [];
for i=1:length(ls)
    set(ls(i),'color', [0, 0, 0, 0.2], 'linewidth', 1);
    x = ls(i).XData; X = [X; x];
    y = ls(i).YData; Y = [Y; y];
end

% Add reference line
x = linspace(0,3.5,10); y = x; 
plot(x,y,'k--', LineWidth=1);

% Add mean data
meanX = mean(X);
meanY = mean(Y);
plot(meanX, meanY, 'Color', config.color_scheme_npg(8,:), 'LineStyle','--','LineWidth',2);

set(gca, ...
'Box'         , 'off'       , ...
'TickDir'     , 'out'       , ...
'TickLength'  , [.01 .01]   , ...
'XColor'      , [.1 .1 .1]  , ...
'YColor'      , [.1 .1 .1]  , ...
'XLim'        , [0, 3.5]     ,...
'YLim'        , [0, 3.5]     ,...
'XTick'       , [0,1,2,3]    ,...
'YTick'       , [0,1,2,3]    ,...
'LineWidth'   , .5          ,...
'XAxisLocation', 'origin'   ,...
'YAxisLocation', 'origin');

xlabel("Intended angle (radians)")
ylabel("Produced agnle (radians)")

exportgraphics(f,config.ResultFolder+"/"+name+"_rgmean.png",'Resolution',300);
exportgraphics(f,config.ResultFolder+"/"+name+"_rgmean.pdf",'Resolution',300, 'ContentType','vector');

%% Plotting the distribution of the single parameters (Fig. S7) 
Parameters = AllHealthyOldResults.estimatedParams;
ParamName = config.ParamName;
for ParamIndx=1:length(ParamName)

    ParamAllConds = []; 

    for TRIAL_FILTER=1:3
        Param = Parameters{TRIAL_FILTER}(:,ParamIndx);
        ParamAllConds = [ParamAllConds,Param];
    end

    ParamAllConds = removeNanRows(ParamAllConds);
    ParamMean = mean(ParamAllConds, 2);
    
    f = figure('visible','off','Position', [100 100 500 500]);

    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)     

    color = config.color_scheme_npg(2,:);

    whisker_value               =   1.5;
    box_lineWidth               =   0.3;
    box_widths_value            =   0.3;
    box_color_transparency      =   0.5; 
    median_lineWidth            =   2;
    median_color                =   'k';
    scatter_jitter_value        =   0.2;
    scatter_markerSize          =   50;
    scatter_marker_edgeColor    =   'k';
    scatter_marker_edgeWidth    =   0.5;
    scatter_color_transparency  =   0.7;         

    bp1 = boxplot(ParamMean, ...
                'Whisker',whisker_value, ...
                'symbol','', ...
                'Color','k', ...
                'Notch','on', ...
                'widths',box_widths_value,...
                'positions', 1);
    set(bp1,'linewidth',box_lineWidth);

    %% Boxplot visual changes
    h = findobj(gca,'Tag','Box'); 
    patch(get(h(1),'XData'),get(h(1),'YData'),color,'FaceAlpha',box_color_transparency);

    %% Median visual changes
    h=findobj(gca,'tag','Median');
    for i = 1:length(h)
        h(i).LineWidth = median_lineWidth;
        h(i).Color = median_color;
    end

    %% Scatter plot for data and mean (MCI positive)
    num_points = length(ParamMean);
    hold on
    x = ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
    scatter(x, ParamMean, scatter_markerSize, ...
            'filled', ...
            'o', ...
            'MarkerEdgeColor',scatter_marker_edgeColor, ...
            'MarkerFaceColor',color, ...
            'MarkerFaceAlpha',scatter_color_transparency,...
            'LineWidth',scatter_marker_edgeWidth); 

    %add errorbar
    mean_ = mean(ParamMean);
    sem_= std(ParamMean)./sqrt(length(ParamMean));
    errorbar(1,mean_,sem_,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18);    
    hold on
    %add mean point
    scatter(1, mean_, 3*scatter_markerSize, 'd',...
            'filled','MarkerEdgeColor','k', ...
            'MarkerFaceColor','w', ...
            'LineWidth',scatter_marker_edgeWidth);

    %% Figure post-processing
    % calculate the Y limits
    alldata = ParamMean;
    maxdata = max(alldata,[],'all');
    mindata = min(alldata, [], 'all');
    
    if ParamName(ParamIndx)=="g2" | ParamName(ParamIndx)=="g3" | ParamName(ParamIndx)=="nu" 
        lowupYlim = [0, 2];
        yticks = [0,1,2,3];
    elseif ParamName(ParamIndx)=="beta"
        lowupYlim = [-0.1, 0.3];
        yticks = [-0.1, 0, 0.1, 0.2, 0.3];
    elseif ParamName(ParamIndx)=="k"
        lowupYlim = [0, 2];
        yticks = [0,0.5,1,1.5,2.0];  
    elseif ParamName(ParamIndx)=="sigma"
        lowupYlim = [0, 1];    
        yticks = [0, 0.5, 1.0];
    end

    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.0 .0 .0], ...
        'YColor'      , [.0 .0 .0], ...
        'XTick'       , (1),... 
        'XLim'        , [0.5, 1.5],...
        'YLim'        , lowupYlim,...   
        'YTick'       , yticks,...
        'XTickLabel'  , {ParamName(ParamIndx)},...
        'LineWidth'   , 1.0        );
    ylabel("Parameter value");

    xL=xlim;
    yL=ylim;

    %Adding reference line and reporting t-tests
    if ParamName(ParamIndx)=="beta"
        [h,p,~,stat]=ttest(ParamMean,0,"Tail","right");
        yline(0,Color='r',LineStyle='--',LineWidth=2);
        text(1.05*xL(1),1.05*yL(2),"t("+num2str(stat.df)+")="+num2str(round(stat.tstat,2))+", p="+num2str(p), 'FontSize', 15)
    elseif ParamName(ParamIndx)=="k"
        [h,p,~,stat]=ttest(ParamMean,1,"Tail","right");
        yline(1,Color='r',LineStyle='--',LineWidth=2); 
        text(1.05*xL(1),1.05*yL(2),"t("+num2str(stat.df)+")="+num2str(round(stat.tstat,2))+", p="+num2str(p), 'FontSize', 15)
    elseif ParamName(ParamIndx)=="g2"
        [h,p,~,stat]=ttest(ParamMean,1);
        yline(1,Color='r',LineStyle='--',LineWidth=2);
        text(1.05*xL(1),1.05*yL(2),"t("+num2str(stat.df)+")="+num2str(round(stat.tstat,2))+", p="+num2str(p), 'FontSize', 15)
    elseif ParamName(ParamIndx)=="g3"
        [h,p,~,stat]=ttest(ParamMean,1,"Tail","left");
        yline(1,Color='r',LineStyle='--',LineWidth=2);
        text(1.05*xL(1),1.05*yL(2),"t("+num2str(stat.df)+")="+num2str(round(stat.tstat,2))+", p="+num2str(p), 'FontSize', 15)
    end

    exportgraphics(f,config.ResultFolder+"/Box_"+ParamName(ParamIndx)+".png",'Resolution',300);
    exportgraphics(f,config.ResultFolder+"/Box_"+ParamName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

end

% Final cleanup to leave workspace as the end of the Preprocessing stage.
% Remove if you want to take a look at the output data.
clearvars -except config YoungControls HealthyControls MCINeg MCIPos MCIUnk

%% ---------------------------------------------------------------------- 
function [Physical_Pos, Mental_Pos, L1_Diff, L2_Diff, L1_Ratio, L2_Ratio, T2, T2_prime, Alpha, T3_prime] = extractAllFinalPoints(GroupResults, config)
    % Function helper to extract all relevant data from the output
    % structure
    Physical_Pos = [];
    Mental_Pos = [];
    L1_Diff = [];
    L2_Diff = [];
    L1_Ratio = [];
    L2_Ratio = [];
    T2 = [];
    T2_prime = [];
    Alpha = {};
    T3_prime = {};
    numID = length(GroupResults.DX{1});
    for ID =1:numID
        L1DiffID = [];
        L2DiffID = [];
        L1RatioID = [];
        L2RatioID = [];
        T2ID = [];
        T2_primeID = [];
        AlphaID = [];
        T3_primeID = [];    
        for cond=1:3
            if ~isempty(GroupResults.DX{cond}{ID}) && ~isnan(GroupResults.estimatedParams{cond}(ID,1))
                numTrails = length(GroupResults.DX{cond}{ID});
                for trial=1:numTrails
                    [phy_p3,men_p3,l1_diff, l2_diff, l1_ratio, l2_ratio, t2, t2_prime, alpha, t3_prime] = VisualizeMenPhyTraj(GroupResults, ID, cond, trial, false, config);

                    Physical_Pos = [Physical_Pos; phy_p3];
                    Mental_Pos = [Mental_Pos; men_p3];
                    L1DiffID = [L1DiffID; l1_diff];
                    L2DiffID = [L2DiffID; l2_diff];

                    L1RatioID = [L1RatioID; l1_ratio];
                    L2RatioID = [L2RatioID; l2_ratio];

                    T2ID = [T2ID; t2];
                    T2_primeID = [T2_primeID; t2_prime];
                    AlphaID = [AlphaID; alpha];
                    T3_primeID = [T3_primeID; t3_prime];
                end
            end
        end

        L1_Diff = [L1_Diff; mean(L1DiffID)];
        L2_Diff = [L2_Diff; mean(L2DiffID)];

        L1_Ratio = [L1_Ratio; mean(L1RatioID)];
        L2_Ratio = [L2_Ratio; mean(L2RatioID)];

        T2 = [T2; mean(T2ID)];
        T2_prime = [T2_prime; mean(T2_primeID)];

        Alpha{ID} = AlphaID;
        T3_prime{ID} = T3_primeID;
    end
end

%% ----------------------------------------------------------------------
function [phy_p3,men_p3_share,  l1_diff, l2_diff, l1_ratio, l2_ratio, theta2, theta2_prime, alpha, theta3_prime]=VisualizeMenPhyTraj(GroupResults, ID, Cond, TrialIdx, doplot, config)
% Visualize mental physical trajectory
    %extract the parameters
    parameters = GroupResults.estimatedParams{Cond}(ID,:);
    cell_params = num2cell(parameters);
    [beta, k, g2, g3, ~, ~] = deal(cell_params{:});

    if isnan(beta)
        error("beta cannot be nan for this plot!")
    end 

    % extract X
    X = GroupResults.X{Cond}{ID}{TrialIdx}; %X here is sufficient to plot the physical trajectory 
    phy_p3 = X(4,:);

    % extract DX, i.e., leg length
    DX              =       GroupResults.DX{Cond}{ID}{TrialIdx}; 
    l1              =       DX(1);
    l2              =       DX(2);

    % extract theta
    Theta           =       GroupResults.THETADX{Cond}{ID}{TrialIdx};
    theta2          =       Theta(2); 
    theta3          =       Theta(3); 

    % find the correct mean return angle based on all trials 
    sampleSize  =   length(GroupResults.DX{1}{ID})+length(GroupResults.DX{2}{ID})+length(GroupResults.DX{3}{ID});
    Alphas      =   zeros(sampleSize,1);
    index       =   0;
    for condddd = 1:3
        trial_num = length(GroupResults.DX{condddd}{ID});
        for trial_id = 1:trial_num
            index = index+1;
            lll_1 = GroupResults.DX{condddd}{ID}{trial_id}(1);
            lll_2 = GroupResults.DX{condddd}{ID}{trial_id}(2);
            thetaaa_2 = GroupResults.THETADX{condddd}{ID}{trial_id}(2);
            
            %calculate the correct return angle
            phy_p1  = [lll_1,0];
            phy_p2  = [lll_1+lll_2*cos(thetaaa_2),lll_2*sin(thetaaa_2)];
            vec1    = phy_p2-phy_p1; vec2 = [0,0]-phy_p2;
            alpha   = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
            alpha   = deg2rad(alpha);%transfer from degree to radians
            alpha   = mod(alpha, 2*pi);  %wrap to (0,2pi)  
            Alphas(index) = alpha;           
        end
    end
    mean_angle = mean(Alphas);

    % extract duration
    durationL1      =       GroupResults.L1Dur{Cond}{ID}{TrialIdx};
    durationL2      =       GroupResults.L2Dur{Cond}{ID}{TrialIdx};

    % calculate the mental trajectory
    men_length1     =       l1*k*(1-exp(-beta*durationL1))/(beta*durationL1)*exp(-beta*durationL2);
    men_p1          =       [men_length1,0];
    
    theta2_prime    =       g2*theta2;

    men_length2     =       l2*k*(1-exp(-beta*durationL2))/(beta*durationL2);
    men_p2          =       [men_length1+men_length2*cos(theta2_prime),men_length2*sin(theta2_prime)];

    %calculate length of mental vector 3
    h               =       norm(men_p2);
    %calculate turn angle of mental vector 3
    vec1            =       men_p2-men_p1; 
    vec2            =       [0,0]-men_p2;
    alpha           =       atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
    alpha           =       deg2rad(alpha);   %transfer from degree to radians
    
    %mental return angle
    sign_alpha      =       sign(alpha);
    theta3_prime    =       g3*abs(alpha)+mean_angle*(1-g3); %reress to mean correct return angle
    theta3_prime    =       sign_alpha*theta3_prime;
        
    
    x3 = h*cos(theta2_prime)*cos(theta3_prime) - h*sin(theta2_prime)*sin(theta3_prime);
    y3 = h*cos(theta2_prime)*sin(theta3_prime) + h*sin(theta2_prime)*cos(theta3_prime);
    men_p3          =       [men_p2(1)+x3, men_p2(2)+y3];
    
    x3_share = h*cos(theta2)*cos(theta3_prime) - h*sin(theta2)*sin(theta3_prime);
    y3_share = h*cos(theta2)*sin(theta3_prime) + h*sin(theta2)*cos(theta3_prime);
    men_p3_share    =       [X(3,1)+x3_share, X(3,2)+y3_share];

    %calculate normalized distance difference
    l1_diff = abs(men_length1-l1)/l1;
    l2_diff = abs(men_length2-l2)/l2;

    %calculate distance ratio
    l1_ratio = men_length1/l1;
    l2_ratio = men_length2/l2;    

    if doplot
        % set figure info
        f = figure('visible','off','Position', [100 100 400 400]);
        %%% Font type and size setting %%%
        % Using Arial as default because all journals normally require the font to
        % be either Arial or Helvetica
        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12) 

        MenColor = config.color_scheme_npg(1,:);
        PhyColor = config.color_scheme_npg(2,:);
    
        %plot mental trajectory
        mentalxy    =       [[0,0];men_p1;men_p2;men_p3];
        m_x         =       mentalxy(:,1); 
        m_y         =       mentalxy(:,2)-0.02;
        plot(m_x', m_y', '--o', 'Markersize', 6,  LineWidth=2, Color=[MenColor,0.8]);
        
        hold on
        %plot physical trajectory
        physicalxy  =       X;
        p_x         =       physicalxy(:,1);
        p_y         =       physicalxy(:,2)+0.02;

        hold on
        plot(p_x', p_y', '.-', 'Markersize', 30,  LineWidth=2, Color=[PhyColor,0.8]);

        set(gca, ...
        'Box'          , 'off'      , ...
        'TickDir'      , 'out'      , ...
        'TickLength'   , [.01 .01]  , ...
        'XColor'       , [.1 .1 .1] , ...
        'YColor'       , [.1 .1 .1] , ...
        'XLim'         , [-2.0, 4.5],...
        'YLim'         , [-2.0, 4.5],... 
        'XTick'        , [-4,0,4]   ,...
        'YTick'        , [-4,0,4]   ,...
        'LineWidth'    , .5         ,...
        'XAxisLocation', 'origin'   ,...
        'YAxisLocation', 'origin');

        title("ID: "+num2str(ID)+" Trial: "+num2str(TrialIdx))

        %% export figure
        exportgraphics(f,config.ResultFolder+"/"+num2str(ID)+"_"+num2str(Cond)+"_"+num2str(TrialIdx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/"+num2str(ID)+"_"+num2str(Cond)+"_"+num2str(TrialIdx)+".pdf",'Resolution',300, 'ContentType','vector');
    end
end