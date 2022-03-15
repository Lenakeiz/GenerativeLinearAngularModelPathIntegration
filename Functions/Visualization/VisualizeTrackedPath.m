%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default');

disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');

%% Plotting some example data
ColorPattern;

Group = MCIPos;
pId = 1;
trialId = 5;
tOffAfterConeReach = 0;

close all;
figure;

% be either Arial or Helvetica
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12)

subplot("Position", [0.05 0.55 0.4 0.4]);

Cone_pos    = Group.FlagPos{1,pId}{trialId,1};
Trig_pos    = Group.TrigPos{1,pId}{trialId,1};
Tracked_pos = Group.Path{1,pId}{trialId,1};
Tracked_pos = array2table(Tracked_pos,"VariableNames",{'Time' 'Pos_X' 'Pos_Y' 'Pos_Z' 'Euler_X' 'Euler_Y' 'Euler_Z'});

hold on;

plCones = plot(Cone_pos(:,1), Cone_pos(:,3),...
    Color=[config.color_scheme_npg(4,:) 1.0],...
    LineStyle="-.",...
    LineWidth=1.3);

plReturn = plot([Cone_pos(3,1);Trig_pos(1,1)], [Cone_pos(3,3);Trig_pos(1,3)],...
    Color=[config.color_scheme_npg(5,:) 1.0],...
    LineStyle=":",...
    LineWidth=1.3);

plPath = plot(Tracked_pos.Pos_X, Tracked_pos.Pos_Z,...
    Color=[config.color_scheme_npg(3,:) 1.0],...
    LineStyle="--",...
    LineWidth=1.3);

hold off;

xlabel('x (m)');
ylabel('y (m)');

legend([plCones plReturn plPath], {'Outbound Path', 'Inbound Path', 'Tracked Data'}, Location="best");

axis square

subplot("Position", [0.55 0.55 0.4 0.4]);
%subplot(2,2,2);

% Calculate the mean time between making cone 1 disappear and spawning cone
% two (same between 2 and 3)

DiffTrigSpawn12 = cell2mat(Group.FlagSpawnTime{1,pId}(:,2)) - cell2mat(Group.FlagTrigTimes{1,pId}(:,1));
DiffTrigSpawn23 = cell2mat(Group.FlagSpawnTime{1,pId}(:,3)) - cell2mat(Group.FlagTrigTimes{1,pId}(:,2));

% Filtering positions
% Excluding tracking before reaching cone 1
Tracked_pos    = Tracked_pos(Tracked_pos.Time > Group.FlagTrigTimes{1,pId}{trialId,1},:);

% Extracting L1 walking path (between flag spawn time 2 and reaching cone 2)
Tracked_pos_L1 = Tracked_pos(Tracked_pos.Time > Group.FlagSpawnTime{1,pId}{trialId,2} &...
                             Tracked_pos.Time < Group.FlagTrigTimes{1,pId}{trialId,2} + mean(DiffTrigSpawn12) + tOffAfterConeReach, :);

% Extracting L1 walking path (between flag spawn time 3 and reaching cone 3)
Tracked_pos_L2 = Tracked_pos(Tracked_pos.Time > Group.FlagSpawnTime{1,pId}{trialId,3} &...
                             Tracked_pos.Time < Group.FlagTrigTimes{1,pId}{trialId,3} + mean(DiffTrigSpawn23) + tOffAfterConeReach, :);


clear tempTrVelL1

hold on;

plCones = plot(Cone_pos(:,1), Cone_pos(:,3),...
    Color=[config.color_scheme_npg(4,:) 1],...
    LineStyle="-.",...
    LineWidth=1.3);

plReturn = plot([Cone_pos(3,1);Trig_pos(1,1)], [Cone_pos(3,3);Trig_pos(1,3)],...
    Color=[config.color_scheme_npg(5,:) 1],...
    LineStyle=":",...
    LineWidth=1.3);

plPathL1 = plot(Tracked_pos_L1.Pos_X, Tracked_pos_L1.Pos_Z,...
    Color=[config.color_scheme_npg(8,:) 1],...
    LineStyle="--",...
    LineWidth=1.3);

plPathL2 = plot(Tracked_pos_L2.Pos_X, Tracked_pos_L2.Pos_Z,...
    Color=[config.color_scheme_npg(9,:) 1],...
    LineStyle="--",...
    LineWidth=1.3);

hold off;

xlabel('x (m)');
ylabel('y (m)');

legend([plCones plReturn plPathL1 plPathL2], {'Outbound Path', 'Inbound Path', 'Tracked Data L1', 'Tracked Data L2'}, Location="best");

axis square

% Calculating ideal direction
s_dir_L1 = (Cone_pos(2,[1 3]) - Cone_pos(1,[1 3])) / norm(Cone_pos(2,[1 3]) - Cone_pos(1,[1 3])); 
s_dir_L2 = (Cone_pos(3,[1 3]) - Cone_pos(2,[1 3])) / norm(Cone_pos(3,[1 3]) - Cone_pos(2,[1 3])); 

%General properties
alpha = 0.9;
smoothWindow = 10;

% Extracting velocity at each time point for segment L1
l1Size = height(Tracked_pos_L1);
Tracked_vel_L1 = array2table(zeros(l1Size,3),"VariableNames",{'Time', 'Vel_X', 'Vel_Z'});
Tracked_vel_L1.Time = Tracked_pos_L1.Time;

for iL1 = 2:l1Size

    Tracked_vel_L1.Vel_X(iL1) = alpha * Tracked_vel_L1.Vel_X(iL1-1) +...
                                (1 - alpha) * (Tracked_pos_L1.Pos_X(iL1) - Tracked_pos_L1.Pos_X(iL1 - 1)) / (Tracked_pos_L1.Time(iL1) - Tracked_pos_L1.Time(iL1-1)); 

    Tracked_vel_L1.Vel_Z(iL1) = alpha * Tracked_vel_L1.Vel_Z(iL1-1) +...
                                (1 - alpha) * (Tracked_pos_L1.Pos_Z(iL1) - Tracked_pos_L1.Pos_Z(iL1 - 1)) / (Tracked_pos_L1.Time(iL1) - Tracked_pos_L1.Time(iL1-1));  

end

% Calculating projecting over the cone2-1 direction
Tracked_vel_L1.Vel_proj = dot([Tracked_vel_L1.Vel_X Tracked_vel_L1.Vel_Z],repmat(s_dir_L1,l1Size,1),2);
Tracked_vel_L1.Smoothed_Vel_proj = smoothdata(Tracked_vel_L1.Vel_proj,'gaussian',smoothWindow);

% Extracting velocity at each time point for segment L2
l2Size = height(Tracked_pos_L2);
Tracked_vel_L2 = array2table(zeros(l2Size,3),"VariableNames",{'Time', 'Vel_X', 'Vel_Z'});
Tracked_vel_L2.Time = Tracked_pos_L2.Time;

for iL1 = 2:l2Size

    Tracked_vel_L2.Vel_X(iL1) = alpha * Tracked_vel_L2.Vel_X(iL1-1) +...
                                (1 - alpha) * (Tracked_pos_L2.Pos_X(iL1) - Tracked_pos_L2.Pos_X(iL1 - 1)) / (Tracked_pos_L2.Time(iL1) - Tracked_pos_L2.Time(iL1-1)); 

    Tracked_vel_L2.Vel_Z(iL1) = alpha * Tracked_vel_L2.Vel_Z(iL1-1) +...
                                (1 - alpha) * (Tracked_pos_L2.Pos_Z(iL1) - Tracked_pos_L2.Pos_Z(iL1 - 1)) / (Tracked_pos_L2.Time(iL1) - Tracked_pos_L2.Time(iL1-1));

end

% Calculating projecting over the cone3-2 direction
Tracked_vel_L2.Vel_proj = dot(table2array(Tracked_vel_L2(:, [2 3])),repmat(s_dir_L2,l2Size,1),2);
Tracked_vel_L2.Smoothed_Vel_proj = smoothdata(Tracked_vel_L2.Vel_proj,'gaussian',smoothWindow);

subplot("Position", [0.05 0.05 0.4 0.4]);

hold on;

pVel = plot(Tracked_vel_L1.Time, vecnorm(table2array(Tracked_vel_L1(:, [2 3])),2,2),...
    Color=[config.color_scheme_npg(8,:) 1],...
    LineStyle="-",...
    LineWidth=1.3);

xlabel('Time (s)');
ylabel('Speed (m/s)');
title('L1');

hold off;

subplot("Position", [0.55 0.05 0.4 0.4]);

hold on;

pVel = plot(Tracked_vel_L2.Time, vecnorm(table2array(Tracked_vel_L2(:, [2 3])),2,2),...
    Color=[config.color_scheme_npg(9,:) 1],...
    LineStyle="-",...
    LineWidth=1.3);

xlabel('Time (s)');
ylabel('Speed (m/s)');
title('L2');

hold off;

figure;

% be either Arial or Helvetica
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',12)

subplot("Position", [0.05 0.05 0.4 0.8]);

hold on;

pProjVel = plot(Tracked_vel_L1.Time, Tracked_vel_L1.Vel_proj,...
    Color=[config.color_scheme_npg(8,:) 1],...
    LineStyle="-",...
    LineWidth=1.3);

pSmoothedProjVel = plot(Tracked_vel_L1.Time, Tracked_vel_L1.Smoothed_Vel_proj,...
    Color=[config.color_scheme_npg(10,:) 1],...
    LineStyle="-",...
    LineWidth=1.3);

xlabel('Time (s)');
ylabel('$\bar{V} \cdot \hat{L1}$','Interpreter','latex');
title('L1');

legend([pProjVel pSmoothedProjVel], {'Projected on L1', 'Gaussian Filtered'}, Location="northwest");

hold off;

subplot("Position", [0.55 0.05 0.4 0.8]);

hold on;

pProjVel = plot(Tracked_vel_L2.Time, Tracked_vel_L2.Vel_proj,...
    Color=[config.color_scheme_npg(9,:) 1],...
    LineStyle="-",...
    LineWidth=1.3);

pSmoothedProjVel = plot(Tracked_vel_L2.Time, Tracked_vel_L2.Smoothed_Vel_proj,...
    Color=[config.color_scheme_npg(2,:) 1],...
    LineStyle="-",...
    LineWidth=1.3);

xlabel('Time (s)');
ylabel('$\bar{V} \cdot \hat{L2}$','Interpreter','latex');
title('L2');

legend([pProjVel pSmoothedProjVel], {'Projected on L2','Guassian filtered'}, Location="northwest");

hold off;