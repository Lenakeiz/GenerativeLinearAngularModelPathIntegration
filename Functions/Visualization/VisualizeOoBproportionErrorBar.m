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

%% statistics of the OoB trial distribution
disp('%%%%%%%%%%%%%%%%%% Statistics of the OoB trial distribution %%%%%%%%%%%%%%%%%%');

YoungRatioMean = zeros(1,3); YoungRatioSem = zeros(1,3); YoungTrialNum = cell(1,3); YoungOoBNum = cell(1,3);
HealthyRatioMean = zeros(1,3); HealthyRatioSem = zeros(1,3); HealthyTrialNum = cell(1,3); HealthyOoBNum = cell(1,3);
MCIPosRatioMean = zeros(1,3); MCIPosRatioSem = zeros(1,3); MCIPosTrialNum = cell(1,3); MCIPosOoBNum = cell(1,3);
MCINegRatioMean = zeros(1,3); MCINegRatioSem = zeros(1,3); MCINegTrialNum = cell(1,3); MCINegOoBNum = cell(1,3);
MCIUnkRatioMean = zeros(1,3); MCIUnkRatioSem = zeros(1,3); MCIUnkTrialNum = cell(1,3); MCIUnkOoBNum = cell(1,3);
%%plot the angle distribution
% 0 all trials, 1 no change, 2 no distal cues, 3 no optic flow

for TRIAL_FILTER=1:3
    %get the distance and angle
    [YoungX, YoungDX, YoungTHETADX, YoungOoBLen, YoungflagOoB] = getDistAngle(YoungControls, TRIAL_FILTER);
    [TrialYoungOoBRatio, TrialNum, OoBNum] = getOoBRatio(YoungflagOoB);
    YoungTrialNum{TRIAL_FILTER} = TrialNum; YoungOoBNum{TRIAL_FILTER}=OoBNum;
    meanYoung = mean(TrialYoungOoBRatio); 
    semYoung =  std(TrialYoungOoBRatio)./sqrt(length(TrialYoungOoBRatio));
    YoungRatioMean(TRIAL_FILTER) = meanYoung; 
    YoungRatioSem(TRIAL_FILTER)=semYoung;

    [HealthyX, HealthyDX, HealthyTHETADX, HealthyOoBLen, HealthyflagOoB] = getDistAngle(HealthyControls, TRIAL_FILTER);
    [TrialHealthyOoBRatio, TrialNum, OoBNum] = getOoBRatio(HealthyflagOoB);
    HealthyTrialNum{TRIAL_FILTER} = TrialNum; HealthyOoBNum{TRIAL_FILTER}=OoBNum;
    meanHealthy = mean(TrialHealthyOoBRatio); 
    semHealthy = std(TrialHealthyOoBRatio)./sqrt(length(TrialHealthyOoBRatio));
    HealthyRatioMean(TRIAL_FILTER) = meanHealthy;
    HealthyRatioSem(TRIAL_FILTER) = semHealthy;

    [MCIPosX, MCIPosDX, MCIPosTHETADX, MCIPosOoBLen, MCIPosflagOoB] = getDistAngle(MCIPos, TRIAL_FILTER);
    [TrialMCIPosOoBRatio, TrialNum, OoBNum] = getOoBRatio(MCIPosflagOoB);
    MCIPosTrialNum{TRIAL_FILTER} = TrialNum; MCIPosOoBNum{TRIAL_FILTER} = OoBNum;
    meanMCIPos = mean(TrialMCIPosOoBRatio); 
    semMCIPos = std(TrialMCIPosOoBRatio)./sqrt(length(TrialMCIPosOoBRatio));
    MCIPosRatioMean(TRIAL_FILTER) = meanMCIPos;
    MCIPosRatioSem(TRIAL_FILTER) = semMCIPos;
    
    [MCINegX, MCINegDX, MCINegTHETADX, MCINegOoBLen, MCINegflagOoB] = getDistAngle(MCINeg, TRIAL_FILTER);
    [TrialMCINegOoBRatio, TrialNum, OoBNum] = getOoBRatio(MCINegflagOoB);
    MCINegTrialNum{TRIAL_FILTER} = TrialNum;MCINegOoBNum{TRIAL_FILTER} = OoBNum;    
    meanMCINeg = mean(TrialMCINegOoBRatio); 
    semMCINeg = std(TrialMCINegOoBRatio)./sqrt(length(TrialMCINegOoBRatio));
    MCINegRatioMean(TRIAL_FILTER) = meanMCINeg;
    MCINegRatioSem(TRIAL_FILTER) = semMCINeg;

    [MCIUnkX, MCIUnkDX, MCIUnkTHETADX, MCIUnkOoBLen, MCIUnkflagOoB] = getDistAngle(Unknown,TRIAL_FILTER);
    [TrialMCIUnkOoBRatio, TrialNum, OoBNum] = getOoBRatio(MCIUnkflagOoB);
    MCIUnkTrialNum{TRIAL_FILTER} = TrialNum;MCIUnkOoBNum{TRIAL_FILTER} = OoBNum;
    meanMCIUnk= mean(TrialMCIUnkOoBRatio); 
    semMCIUnk = std(TrialMCIUnkOoBRatio)./sqrt(length(TrialMCIUnkOoBRatio));
    MCIUnkRatioMean(TRIAL_FILTER) = meanMCIUnk;
    MCIUnkRatioSem(TRIAL_FILTER) = semMCIUnk;
end

Mean = [YoungRatioMean;HealthyRatioMean;MCIPosRatioMean;MCINegRatioMean;MCIUnkRatioMean];  %group X condition
Sem = [YoungRatioSem;HealthyRatioSem;MCIPosRatioSem;MCINegRatioSem;MCIUnkRatioSem]; %group X condition

f = figure('visible','off','Position', [100 100 1000 500]);

b = bar(Mean, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(Mean);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',Mean,Sem,'k','linestyle','none');

legend(b, {'No change' 'No distal cue' 'No optical flow'}, 'Location','northwest');
set(gca, 'XTickLabel', {'Young' 'HealthyOld' 'MCIPos' 'MCINeg' 'MCIUnk'});    
ylabel('Out-of-boundary trial proportion', 'FontSize',20);
exportgraphics(f,"C:/Users/Zilong/Desktop/path integration model/Andrea's matlab code/" + ...
    "PathIntegrationDataEstimationModelling/Output/DataVisualization/OoBErrorBar.png",'Resolution',300);


%% function to get ratio of OoB trials
function [Prop, TrialNum, OoBNum] =getOoBRatio(flagOoB)
%%count the Out-of_boundary trials in a specified group and trial condition
subjectSize = size(flagOoB,2);
%the propotion of OoB trials for each participant and each trial condition
Prop = zeros(1, subjectSize);
TrialNum = zeros(1, subjectSize);
OoBNum = zeros(1, subjectSize);
for subj=1:subjectSize
    subjflagOoB = flagOoB{subj};
    proportion = sum(subjflagOoB)/length(subjflagOoB);
    Prop(subj) = proportion;  
    TrialNum(subj) = length(subjflagOoB);
    OoBNum(subj) = sum(subjflagOoB);
end
end