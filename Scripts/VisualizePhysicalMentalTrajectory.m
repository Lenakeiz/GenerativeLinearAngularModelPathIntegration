%% Visualize physical and mental trajectory 
% create by Zilong, 01/08/2022

%% Preparing the data
VAM_PrepareBaseConfig;

%% Preprocessing the data
VAM_PreprocessData;

%% Preparing the data and Slecting the Model
rng("default")
% config.ModelName        =   "beta_g2_g3_sigma_nu";
% config.ParamName        =   ["beta", "g2", "g3", "sigma", "nu"];

config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];

config.NumParams        =   length(config.ParamName);

% Run the model
VAM;

%% Model run completed, preparing the data for plotting figures
config.ResultFolder     =   pwd + "/Output/ModelFigures/"+config.ModelName+"/VisualizingMentalPhysicalTrajectory";
%create storing folder for trajectory if not exist
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

%% Getting Information from results:
AllYoungResults         =   YoungControls.Results;
AllHealthyOldResults    =   HealthyControls.Results;
AllMCIPosResults        =   MCIPos.Results;
AllMCINegResults        =   MCINeg.Results;
AllMCIUnkResults        =   MCIUnk.Results;

%%
ColorPattern; 

%% visualize one trials
[phy_p3,men_p3] = VisualizeMenPhyTraj(AllHealthyOldResults, 13, 1, 3, true, config); %true means plot the figure

%% visualize some trials
id = 13;
trailNum=4;
%randomly pick a trial in no change condition
%trials = randsample(8,trailNum);
trials = [2,3,5,7];
%trials = [4,6,8,9];

for tN=1:trailNum
    trial_id = trials(tN);
    [phy_p3,men_p3] = VisualizeMenPhyTraj(AllHealthyOldResults, id, 1, trial_id, true, config); %true means plot the figure
end

%% visualize some trials
id = 28;
trailNum=4;
%randomly pick a trial in no change condition
%trials = randsample(8,trailNum);
%trials = [2,3,5,7];
trials = [4,6,8,9];

for tN=1:trailNum
    trial_id = trials(tN);
    [phy_p3,men_p3] = VisualizeMenPhyTraj(AllHealthyOldResults, id, 1, trial_id, true, config); %true means plot the figure
end

%% extract all data informatino
%loop over IDs
%name = "Young";
name = "HealthyOld";
%name = "MCIMerged";

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
    error("choose correct name!")
end

% do stats on return points
f = figure('visible','on','Position', [100 100 500 500]);
%%% Font type and size setting %%%
% Using Arial as default because all journals normally require the font to
% be either Arial or Helvetica
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12) 

MenColor = config.color_scheme_npg(1,:);
PhyColor = config.color_scheme_npg(2,:);

scatter(Physical_Pos(:,1), Physical_Pos(:,2), 'MarkerEdgeColor',PhyColor, ...
    'MarkerFaceColor',PhyColor,'MarkerEdgeAlpha',0.2, 'MarkerFaceAlpha',0.4);
hold on
scatter(Mental_Pos(:,1), Mental_Pos(:,2), 'MarkerEdgeColor',MenColor, ...
    'MarkerFaceColor', MenColor, 'MarkerEdgeAlpha',0.2, 'MarkerFaceAlpha',0.4);

legend('Actual return points', 'Generated return points');

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
%xlabel('X (meters)');
%ylabel('Y (meters)');

exportgraphics(f,config.ResultFolder+"/"+name+".png",'Resolution',300);
exportgraphics(f,config.ResultFolder+"/"+name+".pdf",'Resolution',300, 'ContentType','vector');

%% do stats on return points scatter plot
f = figure('visible','on','Position', [100 100 500 500]);
%%% Font type and size setting %%%
% Using Arial as default because all journals normally require the font to
% be either Arial or Helvetica
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12) 

MenColor = config.color_scheme_npg(1,:);
PhyColor = config.color_scheme_npg(2,:);

Physical_Dist = sqrt(sum(Physical_Pos.^2,2));
Mental_Dist = sqrt(sum(Mental_Pos.^2,2));

%bar scatter plok
scatter(Physical_Dist, Mental_Dist,50, 'MarkerEdgeColor','k', ...
    'MarkerFaceColor','k', 'MarkerEdgeAlpha',0.1, 'MarkerFaceAlpha',0.1)
sh = scatterhist(Physical_Dist,Mental_Dist, 'Kernel','overlay','Location','NorthEast', 'Color','k', 'Direction','out')
x = get(gca,'children');
set(x,'markerfacecolor',[0.9,0.9,0.9])

%scatterhistogram(Physical_Dist,Mental_Dist,'HistogramDisplayStyle','smooth','LineStyle','-')

hold on
% add a linear fit line across origin using the magic backslah operator in
% matlab to do this! see:
% %https://uk.mathworks.com/matlabcentral/answers/392777-how-can-i-plot-the-basic-fitting-line-through-the-origin-0-0
c = Physical_Dist\Mental_Dist;
xc = linspace(0,6,10);
plot(xc,c*xc,'r--', LineWidth=2);

hold on
%add a reference line y=x
x = linspace(0,6,10); y = x; 
plot(x,y,'k--', LineWidth=2);

xlabel("Actual ditsance to origin (meters)")
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
%xlabel('X (meters)');
%ylabel('Y (meters)');

exportgraphics(f,config.ResultFolder+"/"+name+"_distribution.png",'Resolution',300);
exportgraphics(f,config.ResultFolder+"/"+name+"_distribution.pdf",'Resolution',300, 'ContentType','vector');

%% scatter plot of comparing forgetting amount on leg 1 versus leg 2
noNan_L1Diff = L1Diff(~isnan(L1Diff));
noNan_L2Diff = L2Diff(~isnan(L2Diff));

f = figure('visible','on','Position', [100 100, 500, 400]);
%%% Font type and size setting %%%
% Using Arial as default because all journals normally require the font to
% be either Arial or Helvetica
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12) 

Color = config.color_scheme_npg(2,:);

%bar scatter plok
scatter(noNan_L1Diff, noNan_L2Diff,100, 'MarkerEdgeColor','k', ...
    'MarkerFaceColor','k', 'MarkerEdgeAlpha',1.0, 'MarkerFaceAlpha',0.2)
hold on
%add a linear fit line
dlm = fitlm(noNan_L1Diff,noNan_L2Diff,'Intercept',false);
x=0:0.1:1;
yfit = x*dlm.Coefficients{1,1};
plot(x, yfit, 'r--', LineWidth=2)
% ls = lsline;
% ls.LineWidth = 2; 
hold on
%add a reference line y=x
x = linspace(0,1.0,10); y = x; 
plot(x,y,'k--', LineWidth=2);

xlabel("Normalized error for l_1");
ylabel("Normalized error for l_2")

set(gca, ...
'Box'         , 'off'       , ...
'TickDir'     , 'out'       , ...
'TickLength'  , [.01 .01]   , ...
'XColor'      , [.1 .1 .1]  , ...
'YColor'      , [.1 .1 .1]  , ...
'XLim'        , [0, 1.2] ,...
'YLim'        , [0, 1.2] ,... 
'XTick'       , [0,0.5,1.0]    ,...
'YTick'       , [0,0.5,1.0]    ,...
'LineWidth'   , .5          ,...
'XAxisLocation', 'origin'   ,...
'YAxisLocation', 'origin');

exportgraphics(f,config.ResultFolder+"/"+name+"_legdiff.png",'Resolution',300);
exportgraphics(f,config.ResultFolder+"/"+name+"_legdiff.pdf",'Resolution',300, 'ContentType','vector');

%% bar plot of comparing forgetting amount on leg 1 versus leg 2
noNan_L1Ratio = L1Ratio(~isnan(L1Ratio));
noNan_L2Ratio = L2Ratio(~isnan(L2Ratio));

f = figure('visible','on','Position', [100 100, 500, 400]);
%%% Font type and size setting %%%
% Using Arial as default because all journals normally require the font to
% be either Arial or Helvetica
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12) 

% Color1 = config.color_scheme_npg(1,:);
% Color2 = config.color_scheme_npg(2,:);

Color1 = [0.5,0.5,0.5];
Color2 = [0.5,0.5,0.5];


%bar scatter plot
bar(1, mean(noNan_L1Ratio),'FaceColor', Color1, 'EdgeColor', 'k', 'BarWidth',0.4, 'FaceAlpha', 0.8);
hold on
bar(2, mean(noNan_L2Ratio),'FaceColor', Color2, 'EdgeColor', 'k', 'BarWidth',0.4, 'FaceAlpha', 0.8);

%add horizontal connection lines
for i=1:length(noNan_L1Ratio)
    hold on
    plot([1,2], [noNan_L1Ratio(i), noNan_L2Ratio(i)], 'Color',[0,0,0,0.2], LineWidth=1);
end

%scatter plot
hold on
scatter(ones(length(noNan_L1Ratio),1), noNan_L1Ratio,100, 'MarkerEdgeColor','k', ...
    'MarkerFaceColor',Color1, 'MarkerEdgeAlpha',0.2, 'MarkerFaceAlpha',0.2)
hold on
scatter(2*ones(length(noNan_L2Ratio),1), noNan_L2Ratio,100, 'MarkerEdgeColor','k', ...
    'MarkerFaceColor',Color2, 'MarkerEdgeAlpha',0.2, 'MarkerFaceAlpha',0.2)

%link the mean value
hold on
plot([1,2], [mean(noNan_L1Ratio), mean(noNan_L2Ratio)], 'Color',[0,0,0,1.0], LineWidth=2);

hold on
%add a reference line y=1
yl = yline(1);
yl.Color = 'red';
yl.LineWidth = 2;
yl.LineStyle = '--';

[h,p, stats] = ttest(noNan_L1Ratio, noNan_L2Ratio, "Tail", "left")

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

%% second turn barplot
noNan_T2 = T2(~isnan(T2));
noNan_T2_prime = T2_prime(~isnan(T2_prime));

f = figure('visible','on','Position', [100 100, 500, 400]);
%%% Font type and size setting %%%
% Using Arial as default because all journals normally require the font to
% be either Arial or Helvetica
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12) 

MenColor = config.color_scheme_npg(1,:);
PhyColor = config.color_scheme_npg(2,:);

%bar plot
%bar([1,2], [mean(noNan_T2), mean(noNan_T2_prime)],'FaceColor', Color, 'EdgeColor', 'k', 'BarWidth',0.4, 'FaceAlpha', 0.8);
bar(1, mean(noNan_T2),'FaceColor', PhyColor, 'EdgeColor', 'k', 'BarWidth',0.4, 'FaceAlpha', 0.8);
hold on
bar(2, mean(noNan_T2_prime),'FaceColor', MenColor, 'EdgeColor', 'k', 'BarWidth',0.4, 'FaceAlpha', 0.8);

hold on
%error bar
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
%link the mean value
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

%% second turn
% noNan_T2 = T2(~isnan(T2));
% noNan_T2_prime = T2_prime(~isnan(T2_prime));
% 
% f = figure('visible','on','Position', [100 100, 500, 400]);
% %%% Font type and size setting %%%
% % Using Arial as default because all journals normally require the font to
% % be either Arial or Helvetica
% set(0,'DefaultAxesFontName','Arial')
% set(0,'DefaultTextFontName','Arial')
% set(0,'DefaultAxesFontSize',12)
% set(0,'DefaultTextFontSize',12) 
% 
% Color = config.color_scheme_npg(2,:);
% 
% %error bar
% mean_noNan_T2 = mean(noNan_T2);
% sem_noNan_T2 = std(noNan_T2)./sqrt(length(noNan_T2));
% errorbar(1,mean_noNan_T2,sem_noNan_T2,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18); 
% 
% hold on
% mean_noNan_T2_prime = mean(noNan_T2_prime);
% sem_noNan_T2_prime = std(noNan_T2_prime)./sqrt(length(noNan_T2_prime));
% errorbar(2,mean_noNan_T2_prime,sem_noNan_T2_prime,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18); 
% 
% %scatter 
% num_points = size(noNan_T2,1);
% x = 1*ones(num_points,1);
% scatter(x, noNan_T2, 50, 'filled', 'o', 'MarkerEdgeColor',Color, 'MarkerFaceColor',Color, ...
%         'MarkerFaceAlpha',0.2, 'LineWidth',1); 
% hold on
% num_points = size(noNan_T2_prime,1);
% x = 2*ones(num_points,1);
% scatter(x, noNan_T2_prime, 50, 'filled', 'o', 'MarkerEdgeColor',Color, 'MarkerFaceColor',Color, ...
%         'MarkerFaceAlpha',0.5, 'LineWidth',1);
%   
% 
% hold on
% for i=1:num_points
%     plot([1,2], [noNan_T2(i), noNan_T2_prime(i)], 'Color',[0,0,0,0.1], LineWidth=1);
%     hold on
% end
% 
% plot([1,2], [mean(noNan_T2), mean(noNan_T2_prime)], 'Color',[0,0,0,1.0], LineWidth=2);
% 
% set(gca, ...
% 'Box'         , 'off'       , ...
% 'TickDir'     , 'out'       , ...
% 'TickLength'  , [.01 .01]   , ...
% 'XColor'      , [.1 .1 .1]  , ...
% 'YColor'      , [.1 .1 .1]  , ...
% 'XLim'        , [0.5, 2.5] ,...
% 'XTick'       , [1,2]    ,...
% 'XTickLabel'   , {'Physical', 'Mental'},...
% 'LineWidth'   , .5          ,...
% 'XAxisLocation', 'origin'   ,...
% 'YAxisLocation', 'origin');
% 
% exportgraphics(f,config.ResultFolder+"/"+name+"_turn2.png",'Resolution',300);
% exportgraphics(f,config.ResultFolder+"/"+name+"_turn2.pdf",'Resolution',300, 'ContentType','vector');

%% regression to the mean
f = figure('visible','on','Position', [100 100, 500, 400]);
%%% Font type and size setting %%%
% Using Arial as default because all journals normally require the font to
% be either Arial or Helvetica
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

%add linear fitting lines
ls = lsline;
X = [];
Y = [];
for i=1:length(ls)
    %set(ls(i),'color', [CMs(i,:), 0.2], 'linewidth', 1);
    set(ls(i),'color', [0, 0, 0, 0.2], 'linewidth', 1);
    x = ls(i).XData; X = [X; x];
    y = ls(i).YData; Y = [Y; y];
end

%add reference line
%add a reference line y=x
x = linspace(0,3.5,10); y = x; 
plot(x,y,'k--', LineWidth=1);

%add an averaged line
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

%% plot the distribution of all parameters 

Parameters = AllHealthyOldResults.estimatedParams;
ParamName = config.ParamName;
for ParamIndx=1:length(ParamName)

    ParamAllConds = []; 

    for TRIAL_FILTER=1:3
        % extract data
        Param = Parameters{TRIAL_FILTER}(:,ParamIndx);
        ParamAllConds = [ParamAllConds,Param];
    end

    %remove Nan rows (nan coz of 1, removing participants with short walking length; 2, not enough trials for parameter estimation)
    ParamAllConds = removeNanRows(ParamAllConds);
    ParamMean = mean(ParamAllConds, 2);
    
    f = figure('visible','on','Position', [100 100 500 500]);
    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)     
    %%% Color definition %%%
    color = config.color_scheme_npg(2,:);

    %set params
    whisker_value               =   1.5;
    box_lineWidth               =   0.3;
    box_widths_value            =   0.3;
    box_color_transparency      =   0.5; %faceAlpha
    median_lineWidth            =   2;
    median_color                =   'k';
    scatter_jitter_value        =   0.2;
    scatter_markerSize          =   50;
    scatter_marker_edgeColor    =   'k';
    scatter_marker_edgeWidth    =   0.5;
    scatter_color_transparency  =   0.7; %faceAlpha        

    % boxplot for each column in MCIPOs
    bp1 = boxplot(ParamMean, ...
                'Whisker',whisker_value, ...
                'symbol','', ... %symbol ='' making outlier invisible
                'Color','k', ...
                'Notch','on', ...
                'widths',box_widths_value,...
                'positions', 1);
    set(bp1,'linewidth',box_lineWidth);

    % Coloring each box
    h = findobj(gca,'Tag','Box'); 
    %get the MCI box
    patch(get(h(1),'XData'),get(h(1),'YData'),color,'FaceAlpha',box_color_transparency);

    % Adjusting median
    h=findobj(gca,'tag','Median');
    for i = 1:length(h)
        h(i).LineWidth = median_lineWidth;
        h(i).Color = median_color;
    end

    % add scatter plot and the mean of MCIPos
    num_points = length(ParamMean);
    hold on
    x = ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
    scatter(x, ParamMean, scatter_markerSize, ...
            'filled', ...
            'o', ... %marker shape
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

    % Further post-processing the figure
    %calculate the ylim
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
    %add yline and one sample one side ttest
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

%% 
function [Physical_Pos, Mental_Pos, L1_Diff, L2_Diff, L1_Ratio, L2_Ratio, T2, T2_prime, Alpha, T3_prime] = extractAllFinalPoints(GroupResults, config)
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

%% Visualize mental physical trajectory
function [phy_p3,men_p3_share,  l1_diff, l2_diff, l1_ratio, l2_ratio, theta2, theta2_prime, alpha, theta3_prime]=VisualizeMenPhyTraj(GroupResults, ID, Cond, TrialIdx, doplot, config)

    %extract the parameters
    parameters = GroupResults.estimatedParams{Cond}(ID,:);
    cell_params = num2cell(parameters);
    [beta, k, g2, g3, ~, ~] = deal(cell_params{:});

    if isnan(beta)
        error("NAN!")
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
        f = figure('visible','on','Position', [100 100 400 400]);
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

        %plot mental trajectory based on the real location
        %plot([p_x(3), men_p3_share(1)]', [p_y(3), men_p3_share(2)+0.02]', '.-', 'Markersize', 30,  LineWidth=2, Color=[MenColor,0.8]);

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
        %xlabel('X (meters)');
        %ylabel('Y (meters)');
        title("ID: "+num2str(ID)+" Trial: "+num2str(TrialIdx))

        % save figure
        exportgraphics(f,config.ResultFolder+"/"+num2str(ID)+"_"+num2str(Cond)+"_"+num2str(TrialIdx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/"+num2str(ID)+"_"+num2str(Cond)+"_"+num2str(TrialIdx)+".pdf",'Resolution',300, 'ContentType','vector');
    end
end