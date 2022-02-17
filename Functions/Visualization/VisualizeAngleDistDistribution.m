%% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('Data/AllDataErrors2018_V3.mat');

%% Decide whether to add back those Out-of-Boundary trials
% Remember to change the PLOT_FEEDBACK to 0 inside the function inside otherwise it will pause and
% display the paths
disp('%%%%%%%%%%%%%%% TRANSFORME DATA ... %%%%%%%%%%%%%%%');
%Note that I put "CalculateOoBErrors" into TransformPathsOoB to remove the
%repeat calling of "CalculateOoBErrors" here.
YoungControls        = TransformPathsOoB(YoungControls);
HealthyControls      = TransformPathsOoB(HealthyControls);
MCINeg               = TransformPathsOoB(MCINeg);
MCIPos               = TransformPathsOoB(MCIPos);
Unknown              = TransformPathsOoB(Unknown);

%% plot the distance distribution
disp('%%%%%%%%%%%%%%%%%% Plot distance distribution %%%%%%%%%%%%%%%%%%');
legidx = 1;
%%plot the angle distribution
% 0 all trials, 1 no change, 2 no distal cues, 3 no optic flow
for TRIAL_FILTER=0:3

    % plot histogram of all subjects' third turn in subplots
    f = figure('visible','off','Position', [100 100 600 1000]);

    subplot(5,1,1)
    %get the distance and angle
    [YoungX, YoungDX, YoungTHETADX, YoungOoBLen, YoungflagOoB] = getDistAngle(YoungControls, TRIAL_FILTER);
    [YoungDist, YoungflagOoBs] = get_distance3_from_groupdata(YoungDX, YoungOoBLen, YoungflagOoB, legidx);
    %non OoB trials
    histogram(YoungDist(~YoungflagOoBs), 30, 'FaceColor','m', 'Facealpha',0.5, 'EdgeColor', 'm', 'EdgeAlpha', 0.5);
    hold on
    %OoB trials
    histogram(YoungDist(YoungflagOoBs), 30, 'FaceColor','k', 'Facealpha',0.5, 'EdgeColor', 'k', 'EdgeAlpha', 0.5);
    ylabel('frequency', 'FontSize', 15); xlabel('Young', 'FontSize', 15);
    xticks([ 0, 2, 4, 6]);
    xlim([0, 6]);    

    subplot(5,1,2)
    %get the distance and angle
    [HealthyX, HealthyDX, HealthyTHETADX, HealthyOoBLen, HealthyflagOoB] = getDistAngle(HealthyControls, TRIAL_FILTER);
    [HealthyDist, HealthyflagOoBs] = get_distance3_from_groupdata(HealthyDX, HealthyOoBLen, HealthyflagOoB, legidx);

    if TRIAL_FILTER==0 HealthyABNMatrix=checkabnormal(HealthyDX); end
    
    %non OoB trials
    histogram(HealthyDist(~HealthyflagOoBs), 30, 'FaceColor','m', 'Facealpha',0.5, 'EdgeColor', 'm', 'EdgeAlpha', 0.5);
    hold on
    %OoB trials
    histogram(HealthyDist(HealthyflagOoBs), 30, 'FaceColor','k', 'Facealpha',0.5, 'EdgeColor', 'k', 'EdgeAlpha', 0.5);
    ylabel('frequency', 'FontSize', 15);xlabel('HealthyOld', 'FontSize', 15);
    xticks([ 0, 2, 4, 6]);
    xlim([0, 6]);     

    subplot(5,1,3)
    %get the distance and angle
    [MCIPosX, MCIPosDX, MCIPosTHETADX, MCIPosOoBLen, MCIPosflagOoB] = getDistAngle(MCIPos, TRIAL_FILTER);
    [MCIPosDist, MCIPosflagOoBs] = get_distance3_from_groupdata(MCIPosDX, MCIPosOoBLen, MCIPosflagOoB, legidx);
    %non OoB trials
    histogram(MCIPosDist(~MCIPosflagOoBs), 30, 'FaceColor','m', 'Facealpha',0.5, 'EdgeColor', 'm', 'EdgeAlpha', 0.5);
    hold on
    %OoB trials
    histogram(MCIPosDist(MCIPosflagOoBs), 30, 'FaceColor','k', 'Facealpha',0.5, 'EdgeColor', 'k', 'EdgeAlpha', 0.5);
    ylabel('frequency', 'FontSize', 15);xlabel('MCIPos', 'FontSize', 15);
    xticks([ 0, 2, 4, 6]);
    xlim([0, 6]);   

    subplot(5,1,4)
    %get the distance and angle
    [MCINegX, MCINegDX, MCINegTHETADX, MCINegOoBLen, MCINegflagOoB] = getDistAngle(MCINeg, TRIAL_FILTER);
    [MCINegDist, MCINegflagOoBs] = get_distance3_from_groupdata(MCINegDX, MCINegOoBLen, MCINegflagOoB, legidx);
    %non OoB trials
    histogram(MCINegDist(~MCINegflagOoBs), 30, 'FaceColor','m', 'Facealpha',0.5, 'EdgeColor', 'm', 'EdgeAlpha', 0.5);
    hold on
    %OoB trials
    histogram(MCINegDist(MCINegflagOoBs), 30, 'FaceColor','k', 'Facealpha',0.5, 'EdgeColor', 'k', 'EdgeAlpha', 0.5);
    ylabel('frequency', 'FontSize', 15);xlabel('MCINeg', 'FontSize', 15);
    xticks([ 0, 2, 4, 6]);
    xlim([0, 6]);     

    subplot(5,1,5)
    %get the distance and angle
    [MCIUnkX, MCIUnkDX, MCIUnkTHETADX, MCIUnkOoBLen, MCIUnkflagOoB] = getDistAngle(Unknown,TRIAL_FILTER);
    [MCIUnkDist, MCIUnkflagOoBs] = get_distance3_from_groupdata(MCIUnkDX, MCIUnkOoBLen, MCIUnkflagOoB, legidx);
    if TRIAL_FILTER==0 MCIUnkABNMatrix=checkabnormal(MCIUnkDX); end
    %non OoB trials
    histogram(MCIUnkDist(~MCIUnkflagOoBs), 30, 'FaceColor','m', 'Facealpha',0.5, 'EdgeColor', 'm', 'EdgeAlpha', 0.5);
    hold on
    %OoB trials
    histogram(MCIUnkDist(MCIUnkflagOoBs), 30, 'FaceColor','k', 'Facealpha',0.5, 'EdgeColor', 'k', 'EdgeAlpha', 0.5);
    ylabel('frequency', 'FontSize', 15);xlabel('MCIUnk', 'FontSize', 15);
    xticks([ 0, 2, 4, 6]);
    xlim([0, 6]);    

    sgtitle("Distance Distribution of Leg " +num2str(legidx), 'FontSize', 20);
    
    exportgraphics(f,"C:/Users/Zilong/Desktop/path integration model/Andrea's matlab code/" + ...
    "PathIntegrationDataEstimationModelling/Output/DataVisualization/legidx_"+num2str(legidx)+ ...
    "/Distance3_TrialType_" + num2str(TRIAL_FILTER)+".png",'Resolution',300);

    close(f);
end

%%
f = figure('visible','off','Position', [100 100 1000 1000]);
image(HealthyABNMatrix*255);
[rows, columns] = size(HealthyABNMatrix);
xlabel('trial index', 'FontSize', 40); ylabel('subject index', 'FontSize', 40);
exportgraphics(f,"C:/Users/Zilong/Desktop/path integration model/Andrea's matlab code/" + ...
    "PathIntegrationDataEstimationModelling/Output/DataVisualization/HealthyOldABN.png",'Resolution',300);
close(f);

f = figure('visible','off','Position', [100 100 1000 1000]);
image(MCIUnkABNMatrix*255);
xlabel('trial index', 'FontSize', 40); ylabel('subject index', 'FontSize', 40);
exportgraphics(f,"C:/Users/Zilong/Desktop/path integration model/Andrea's matlab code/" + ...
    "PathIntegrationDataEstimationModelling/Output/DataVisualization/MCIUnkABN.png",'Resolution',300);
close(f);
%% plot the angle distribution
disp('%%%%%%%%%%%%%%%%%% Plot angle distribution %%%%%%%%%%%%%%%%%%');
% 0 all trials, 1 no change, 2 no distal cues, 3 no optic flow
for TRIAL_FILTER=0:3

    % plot histogram of all subjects' third turn in subplots
    f = figure('visible','off','Position', [100 100 600 1000]);

    subplot(5,1,1)
    %get the distance and angle
    [YoungX, YoungDX, YoungTHETADX, YoungOoBLen, YoungflagOoB] = getDistAngle(YoungControls, TRIAL_FILTER);
    [YoungAngle, YoungflagOoBs] = get_angle3_from_groupdata(YoungTHETADX, YoungflagOoB, legidx);
    %non OoB trials
    histogram(YoungAngle(~YoungflagOoBs), 30, 'FaceColor','m', 'Facealpha',0.5, 'EdgeColor', 'm', 'EdgeAlpha', 0.5);
    hold on
    %OoB trials
    histogram(YoungAngle(YoungflagOoBs), 30, 'FaceColor','k', 'Facealpha',0.5, 'EdgeColor', 'k', 'EdgeAlpha', 0.5);
    ylabel('frequency', 'FontSize', 15); xlabel('Young', 'FontSize', 15);
    xticks([ 0, pi/2, pi, 3/2*pi, 2*pi]);
    xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
    xlim([0, 2*pi]);    

    subplot(5,1,2)
    %get the distance and angle
    [HealthyX, HealthyDX, HealthyTHETADX, HealthyOoBLen, HealthyflagOoB] = getDistAngle(HealthyControls, TRIAL_FILTER);
    [HealthyAngle, HealthyflagOoBs] = get_angle3_from_groupdata(HealthyTHETADX, HealthyflagOoB, legidx);
    %non OoB trials
    histogram(HealthyAngle(~HealthyflagOoBs), 30, 'FaceColor','m', 'Facealpha',0.5, 'EdgeColor', 'm', 'EdgeAlpha', 0.5);
    hold on
    %OoB trials
    histogram(HealthyAngle(HealthyflagOoBs), 30, 'FaceColor','k', 'Facealpha',0.5, 'EdgeColor', 'k', 'EdgeAlpha', 0.5);
    ylabel('frequency', 'FontSize', 15);xlabel('HealthyOld', 'FontSize', 15);
    xticks([ 0, pi/2, pi, 3/2*pi, 2*pi]);
    xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
    xlim([0, 2*pi]);    

    subplot(5,1,3)
    %get the distance and angle
    [MCIPosX, MCIPosDX, MCIPosTHETADX, MCIPosOoBLen, MCIPosflagOoB] = getDistAngle(MCIPos, TRIAL_FILTER);
    [MCIPosAngle, MCIPosflagOoBs] = get_angle3_from_groupdata(MCIPosTHETADX, MCIPosflagOoB, legidx);
    %non OoB trials
    histogram(MCIPosAngle(~MCIPosflagOoBs), 30, 'FaceColor','m', 'Facealpha',0.5, 'EdgeColor', 'm', 'EdgeAlpha', 0.5);
    hold on
    %OoB trials
    histogram(MCIPosAngle(MCIPosflagOoBs), 30, 'FaceColor','k', 'Facealpha',0.5, 'EdgeColor', 'k', 'EdgeAlpha', 0.5);
    ylabel('frequency', 'FontSize', 15);xlabel('MCIPos', 'FontSize', 15);
    xticks([ 0, pi/2, pi, 3/2*pi, 2*pi]);
    xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
    xlim([0, 2*pi]);    

    subplot(5,1,4)
    %get the distance and angle
    [MCINegX, MCINegDX, MCINegTHETADX, MCINegOoBLen, MCINegflagOoB] = getDistAngle(MCINeg, TRIAL_FILTER);
    [MCINegAngle, MCINegflagOoBs] = get_angle3_from_groupdata(MCINegTHETADX, MCINegflagOoB, legidx);
    %non OoB trials
    histogram(MCINegAngle(~MCINegflagOoBs), 30, 'FaceColor','m', 'Facealpha',0.5, 'EdgeColor', 'm', 'EdgeAlpha', 0.5);
    hold on
    %OoB trials
    histogram(MCINegAngle(MCINegflagOoBs), 30, 'FaceColor','k', 'Facealpha',0.5, 'EdgeColor', 'k', 'EdgeAlpha', 0.5);
    ylabel('frequency', 'FontSize', 15);xlabel('MCINeg', 'FontSize', 15);
    xticks([ 0, pi/2, pi, 3/2*pi, 2*pi]);
    xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
    xlim([0, 2*pi]);    

    subplot(5,1,5)
    %get the distance and angle
    [MCIUnkX, MCIUnkDX, MCIUnkTHETADX, MCIUnkOoBLen, MCIUnkflagOoB] = getDistAngle(Unknown,TRIAL_FILTER);
    [MCIUnkAngle, MCIUnkflagOoBs] = get_angle3_from_groupdata(MCIUnkTHETADX, MCIUnkflagOoB, legidx);
    %non OoB trials
    histogram(MCIUnkAngle(~MCIUnkflagOoBs), 30, 'FaceColor','m', 'Facealpha',0.5, 'EdgeColor', 'm', 'EdgeAlpha', 0.5);
    hold on
    %OoB trials
    histogram(MCIUnkAngle(MCIUnkflagOoBs), 30, 'FaceColor','k', 'Facealpha',0.5, 'EdgeColor', 'k', 'EdgeAlpha', 0.5);
    ylabel('frequency', 'FontSize', 15);xlabel('MCIUnk', 'FontSize', 15);
    xticks([ 0, pi/2, pi, 3/2*pi, 2*pi]);
    xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
    xlim([0, 2*pi]);    
    
    sgtitle("Angle Distribution of Leg " +num2str(legidx), 'FontSize', 20);

    exportgraphics(f,"C:/Users/Zilong/Desktop/path integration model/Andrea's matlab code/" + ...
    "PathIntegrationDataEstimationModelling/Output/DataVisualization/legidx_"+num2str(legidx)+ ...
    "/Angle3_TrialType_" + num2str(TRIAL_FILTER)+".png",'Resolution',300);

    close(f);
end


%% function to get distance 3
function [all_dists,all_flagOoBs]=get_distance3_from_groupdata(DX, OoBLen, flagOoB, legidx)
%%plot angle 3 distibution
% DX is a cell structure containing the segment of each trial
% flagOoB is the flag of Out-of-Boundary trials, with 0 non-OoB trials 1 those OoB trials

subjectSize = size(DX,2);

all_dists = [];
all_flagOoBs = [];

for subj=1:subjectSize
    subjDX = DX{subj};
    subjOoBLen = OoBLen{subj};
    subjflagOoB = flagOoB{subj};

    sampleSize = size(subjDX,2);
    Dists = zeros(1,sampleSize);
    for tr_id=1:sampleSize
        if subjflagOoB(tr_id)==0
            %non-OoB trial
            dist3 = subjDX{tr_id}(legidx);
        else
            %OoB trial
            if legidx ==3
                %at the 3rd leg, using the hit-boundary distance 
                dist3 = subjOoBLen(tr_id);
            else
                %at the 1st, 2nd leg, using the taken distance
                dist3 = subjDX{tr_id}(legidx);
            end
        end
        Dists(tr_id) = dist3;
    end
    all_dists = [all_dists,Dists];
    all_flagOoBs = [all_flagOoBs, subjflagOoB]; 
    %turn to logcial value for indexing
    all_flagOoBs = logical(all_flagOoBs);    
end
end



%% function to get angle3
function [all_angles,all_flagOoBs]=get_angle3_from_groupdata(THETADX, flagOoB, legidx)
%%plot angle 3 distibution
% THETAX is the turning angle (wrong at the moment)
% flagOoB is the flag of Out-of-Boundary trials, with 0 non-OoB trials 1 those OoB trials

subjectSize = size(THETADX,2);

all_angles = [];
all_flagOoBs = [];

for subj=1:subjectSize
    subjTHETADX = THETADX{subj};
    subjflagOoB = flagOoB{subj};

    sampleSize = size(subjTHETADX,2);
    Angles = zeros(1,sampleSize);
    for tr_id=1:sampleSize
        turn3_angle = subjTHETADX{tr_id}(legidx);
        Angles(tr_id) = turn3_angle;
    end
    all_angles = [all_angles,Angles];
    all_flagOoBs = [all_flagOoBs, subjflagOoB]; 
    %turn to logcial value for indexing
    all_flagOoBs = logical(all_flagOoBs);    
end
end

function ABNMatrix=checkabnormal(DX)
    subjectSize = size(DX,2);
    sampleSize = size(DX{1},2);
    ABNMatrix = zeros(subjectSize, sampleSize);
    for subj = 1:subjectSize
        subjDX = DX{subj};
        for tr_id = 1:sampleSize
            flag = subjDX{tr_id}(1)<2;
            ABNMatrix(subj,tr_id) = flag;
        end
    end
end
