%% FittingData 
% This script will load the data, rotate the paths and call the solver to find the parameters 
% that minimize the mean distance to generated x3' which matchs the physical leg 3

%% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V2.mat');

%% Decide whether to add back those Out-of-Boundary trials
% Remember to change the PLOT_FEEDBACK to 0 inside the function inside otherwise it will pause and
% display the paths
disp('%%%%%%%%%%%%%%% TRANSFORME DATA ... %%%%%%%%%%%%%%%');
%Note that I put "CalculateOoBErrors" into TransformPathsOoB to remove the
%repeat calling of "CalculateOoBErrors" here.
YoungControls        = TransformPaths(YoungControls);
HealthyControls      = TransformPaths(HealthyControls);
MCINeg               = TransformPaths(MCINeg);
MCIPos               = TransformPaths(MCIPos);
Unknown              = TransformPaths(Unknown);

%% plot the distance distribution
disp('%%%%%%%%%%%%%%%%%% Plot distance distribution %%%%%%%%%%%%%%%%%%');
legidx = 2;
%%plot the angle distribution
% 0 all trials, 1 no change, 2 no distal cues, 3 no optic flow
for TRIAL_FILTER=0:3

    % plot histogram of all subjects' third turn in subplots
    f = figure('visible','off','Position', [100 100 1000 800]);

    subplot(5,1,1)
    %get the distance and angle
    [YoungX, YoungDX, YoungTHETADX] = getDistAngle_v2(YoungControls, TRIAL_FILTER);
    [YoungDist] = get_distance3_from_groupdata(YoungDX, legidx);
    histogram(YoungDist, 30, 'FaceColor','m', 'facealpha',0.5);
    ylabel('frequency'); xlabel('Young', 'FontSize', 10);
    xticks([ 0, 2, 4, 6]);
    xlim([0, 6]);    

    subplot(5,1,2)
    %get the distance and angle
    [HealthyX, HealthyDX, HealthyTHETADX] = getDistAngle_v2(HealthyControls, TRIAL_FILTER);
    [HealthyDist] = get_distance3_from_groupdata(HealthyDX, legidx);
    %non OoB trials
    histogram(HealthyDist, 30, 'FaceColor','m', 'facealpha',0.5);
    ylabel('frequency');xlabel('HealthyOld', 'FontSize', 10);
    xticks([ 0, 2, 4, 6]);
    xlim([0, 6]);     

    subplot(5,1,3)
    %get the distance and angle
    [MCIPosX, MCIPosDX, MCIPosTHETADX] = getDistAngle_v2(MCIPos, TRIAL_FILTER);
    [MCIPosDist] = get_distance3_from_groupdata(MCIPosDX, legidx);
    histogram(MCIPosDist, 30, 'FaceColor','m', 'facealpha',0.5);
    ylabel('frequency');xlabel('MCIPos', 'FontSize', 10);
    xticks([ 0, 2, 4, 6]);
    xlim([0, 6]);   

    subplot(5,1,4)
    %get the distance and angle
    [MCINegX, MCINegDX, MCINegTHETADX] = getDistAngle_v2(MCINeg, TRIAL_FILTER);
    [MCINegDist] = get_distance3_from_groupdata(MCINegDX, legidx);
    histogram(MCINegDist, 30, 'FaceColor','m', 'facealpha',0.5);
    ylabel('frequency');xlabel('MCINeg', 'FontSize', 10);
    xticks([ 0, 2, 4, 6]);
    xlim([0, 6]);     

    subplot(5,1,5)
    %get the distance and angle
    [MCIUnkX, MCIUnkDX, MCIUnkTHETADX] = getDistAngle_v2(Unknown,TRIAL_FILTER);
    [MCIUnkDist] = get_distance3_from_groupdata(MCIUnkDX, legidx);
    histogram(MCIUnkDist, 30, 'FaceColor','m', 'facealpha',0.5);
    ylabel('frequency');xlabel('MCIUnk', 'FontSize', 10);
    xticks([ 0, 2, 4, 6]);
    xlim([0, 6]);    

    sgtitle("Leg " +num2str(legidx)+" distribution of all participants");
    
    exportgraphics(f,"C:/Users/Zilong/Desktop/path integration model/Andrea's matlab code/" + ...
    "PathIntegrationDataEstimationModelling/Output/DataVisualization_v2/legidx_"+num2str(legidx)+ ...
    "/Distance3_TrialType_" + num2str(TRIAL_FILTER)+".png",'Resolution',300);

    close(f);
end


%% plot the angle distribution
disp('%%%%%%%%%%%%%%%%%% Plot angle distribution %%%%%%%%%%%%%%%%%%');
% 0 all trials, 1 no change, 2 no distal cues, 3 no optic flow
for TRIAL_FILTER=0:3

    % plot histogram of all subjects' third turn in subplots
    f = figure('visible','off','Position', [100 100 1000 800]);

    subplot(5,1,1)
    %get the distance and angle
    [YoungX, YoungDX, YoungTHETADX] = getDistAngle_v2(YoungControls, TRIAL_FILTER);
    [YoungAngle] = get_angle3_from_groupdata(YoungTHETADX, legidx);
    histogram(YoungAngle, 30, 'FaceColor','m', 'facealpha',0.5);
    ylabel('frequency'); xlabel('Young', 'FontSize', 10);
    xticks([ 0, pi/2, pi, 3/2*pi, 2*pi]);
    xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
    xlim([0, 2*pi]);    

    subplot(5,1,2)
    %get the distance and angle
    [HealthyX, HealthyDX, HealthyTHETADX] = getDistAngle_v2(HealthyControls, TRIAL_FILTER);
    [HealthyAngle] = get_angle3_from_groupdata(HealthyTHETADX, legidx);
    histogram(HealthyAngle, 30, 'FaceColor','m', 'facealpha',0.5);
    ylabel('frequency');xlabel('HealthyOld', 'FontSize', 10);
    xticks([ 0, pi/2, pi, 3/2*pi, 2*pi]);
    xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
    xlim([0, 2*pi]);    

    subplot(5,1,3)
    %get the distance and angle
    [MCIPosX, MCIPosDX, MCIPosTHETADX] = getDistAngle_v2(MCIPos, TRIAL_FILTER);
    [MCIPosAngle] = get_angle3_from_groupdata(MCIPosTHETADX, legidx);
    histogram(MCIPosAngle, 30, 'FaceColor','m', 'facealpha',0.5);
    ylabel('frequency');xlabel('MCIPos', 'FontSize', 10);
    xticks([ 0, pi/2, pi, 3/2*pi, 2*pi]);
    xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
    xlim([0, 2*pi]);    

    subplot(5,1,4)
    %get the distance and angle
    [MCINegX, MCINegDX, MCINegTHETADX] = getDistAngle_v2(MCINeg, TRIAL_FILTER);
    [MCINegAngle] = get_angle3_from_groupdata(MCINegTHETADX, legidx);
    histogram(MCINegAngle, 30, 'FaceColor','m', 'facealpha',0.5);
    ylabel('frequency');xlabel('MCINeg', 'FontSize', 10);
    xticks([ 0, pi/2, pi, 3/2*pi, 2*pi]);
    xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
    xlim([0, 2*pi]);    

    subplot(5,1,5)
    %get the distance and angle
    [MCIUnkX, MCIUnkDX, MCIUnkTHETADX] = getDistAngle_v2(Unknown,TRIAL_FILTER);
    [MCIUnkAngle] = get_angle3_from_groupdata(MCIUnkTHETADX, legidx);
    histogram(MCIUnkAngle, 30, 'FaceColor','m', 'facealpha',0.5);
    ylabel('frequency');xlabel('MCIUnk', 'FontSize', 10);
    xticks([ 0, pi/2, pi, 3/2*pi, 2*pi]);
    xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
    xlim([0, 2*pi]);    
    
    sgtitle("Angle " +num2str(legidx)+" distribution of all participants")

    exportgraphics(f,"C:/Users/Zilong/Desktop/path integration model/Andrea's matlab code/" + ...
    "PathIntegrationDataEstimationModelling/Output/DataVisualization_v2/legidx_"+num2str(legidx)+ ...
    "/Angle3_TrialType_" + num2str(TRIAL_FILTER)+".png",'Resolution',300);

    close(f);
end


%% function to get distance 3
function [all_dists]=get_distance3_from_groupdata(DX, legidx)
%%plot angle 3 distibution
% DX is a cell structure containing the segment of each trial

subjectSize = size(DX,2);

all_dists = [];

for subj=1:subjectSize
    subjDX = DX{subj};

    sampleSize = size(subjDX,2);
    Dists = zeros(1,sampleSize);
    for tr_id=1:sampleSize
        Dists(tr_id) = subjDX{tr_id}(legidx);
    end
    all_dists = [all_dists,Dists];
end
end



%% function to get angle3
function [all_angles]=get_angle3_from_groupdata(THETADX, legidx)
%%plot angle 3 distibution
% THETAX is the turning angle (wrong at the moment)

subjectSize = size(THETADX,2);

all_angles = [];

for subj=1:subjectSize
    subjTHETADX = THETADX{subj};

    sampleSize = size(subjTHETADX,2);
    Angles = zeros(1,sampleSize);
    for tr_id=1:sampleSize
        turn3_angle = subjTHETADX{tr_id}(legidx);
        Angles(tr_id) = turn3_angle;
    end
    all_angles = [all_angles,Angles];
end
end
