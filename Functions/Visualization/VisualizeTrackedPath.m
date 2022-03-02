%% Cleaning variables
clearvars; clear all; close all; clc;
rng('default');

disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');

%% Plotting some example data
ColorPattern;

Group = HealthyControls;
pId = 15;
trialId = 3;

close all;
figure;

subplot("Position", [0.05 0.05 0.4 0.95]);

Cone_pos    = Group.FlagPos{1,pId}{trialId,1};
Trig_pos    = Group.TrigPos{1,pId}{trialId,1};
Tracked_pos = Group.Path{1,pId}{trialId,1};
Tracked_pos = array2table(Tracked_pos,"VariableNames",{'Time' 'Pos_X' 'Pos_Y' 'Pos_Z' 'Euler_X' 'Euler_Y' 'Euler_Z'});

hold on;

plCones = plot(Cone_pos(:,1), Cone_pos(:,3),...
    Color=[config.color_scheme_npg(2,:) 1.0],...
    LineStyle="-.",...
    LineWidth=1.3);

plReturn = plot([Cone_pos(3,1);Trig_pos(1,1)], [Cone_pos(3,3);Trig_pos(1,3)],...
    Color=[config.color_scheme_npg(1,:) 1.0],...
    LineStyle=":",...
    LineWidth=1.3);

plPath = plot(Tracked_pos.Pos_X, Tracked_pos.Pos_Z,...
    Color=[config.color_scheme_npg(3,:) 1.0],...
    LineStyle="--",...
    LineWidth=1.3);

hold off;

xlabel('x (m)');
ylabel('y (m)');

legend([plCones plReturn plPath], {'Outbound Path', 'Inbound Path', 'Tracked Data'}, Location="northoutside");

axis square

subplot("Position", [0.55 0.05 0.4 0.95]);

% Calculate the mean time between making cone 1 disappear and spawning cone
% two (same between 2 and 3)

DiffTrigSpawn12 = cell2mat(Group.FlagSpawnTime{1,pId}(:,2)) - cell2mat(Group.FlagTrigTimes{1,pId}(:,1));
DiffTrigSpawn23 = cell2mat(Group.FlagSpawnTime{1,pId}(:,3)) - cell2mat(Group.FlagTrigTimes{1,pId}(:,2));

% Filtering positions
% Excluding tracking before reaching cone 1
Tracked_pos    = Tracked_pos(Tracked_pos.Time > Group.FlagTrigTimes{1,pId}{trialId,1},:);

% Extracting L1 walking path (between flag spawn time 2 and reaching cone 2)
Tracked_pos_L1 = Tracked_pos(Tracked_pos.Time > Group.FlagSpawnTime{1,pId}{trialId,2} &...
                             Tracked_pos.Time < Group.FlagTrigTimes{1,pId}{trialId,2} + mean(DiffTrigSpawn12)/2 , :);

% Extracting L1 walking path (between flag spawn time 2 and reaching cone 2)
Tracked_pos_L2 = Tracked_pos(Tracked_pos.Time > Group.FlagSpawnTime{1,pId}{trialId,3} &...
                             Tracked_pos.Time < Group.FlagTrigTimes{1,pId}{trialId,3} + mean(DiffTrigSpawn23)/2, :);

% Extracting velocity at each time  point

hold on;

plCones = plot(Cone_pos(:,1), Cone_pos(:,3),...
    Color=[config.color_scheme_npg(2,:) 1],...
    LineStyle="-.",...
    LineWidth=1.3);

plReturn = plot([Cone_pos(3,1);Trig_pos(1,1)], [Cone_pos(3,3);Trig_pos(1,3)],...
    Color=[config.color_scheme_npg(1,:) 1],...
    LineStyle=":",...
    LineWidth=1.3);

plPathL1 = plot(Tracked_pos_L1.Pos_X, Tracked_pos_L1.Pos_Z,...
    Color=[config.color_scheme_npg(3,:) 1],...
    LineStyle="--",...
    LineWidth=1.3);

plPathL2 = plot(Tracked_pos_L2.Pos_X, Tracked_pos_L2.Pos_Z,...
    Color=[config.color_scheme_npg(4,:) 1],...
    LineStyle="--",...
    LineWidth=1.3);

hold off;

xlabel('x (m)');
ylabel('y (m)');

legend([plCones plReturn plPathL1 plPathL2], {'Outbound Path', 'Inbound Path', 'Tracked Data L1', 'Tracked Data L2'}, Location="northoutside");

axis square

clear Cone_pos Trig_pos Tracked_pos plCones plReturn plPath pId trialId
