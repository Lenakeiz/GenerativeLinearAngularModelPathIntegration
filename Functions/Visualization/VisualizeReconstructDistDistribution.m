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
Data = cell(1,5);
DataName = ["Young", "HealthyOld", "MCIPos", "MCINeg", "MCIUnk"];
YoungControls        = TransformPathsOoB(YoungControls); Data{1,1} = YoungControls;
HealthyControls      = TransformPathsOoB(HealthyControls);Data{1,2} = HealthyControls;
MCINeg               = TransformPathsOoB(MCINeg);Data{1,3} = MCINeg;
MCIPos               = TransformPathsOoB(MCIPos);Data{1,4} = MCIPos;
Unknown              = TransformPathsOoB(Unknown);Data{1,5} = Unknown;

%% plot the distance distribution
disp('%%%%%%%%%%%%%%%%%% Plot distance distribution %%%%%%%%%%%%%%%%%%');
%%plot the angle distribution
% 0 all trials, 1 no change, 2 no distal cues, 3 no optic flow
trialname = ["No change", "No distal cue", "No optical flow"];

for groupid = 1:5
    dat = Data{groupid};
    datname = DataName(groupid);

    % plot histogram of all subjects' third turn in subplots
    f = figure('visible','off','Position', [100 100 1000 300]);

    for TRIAL_FILTER=1:3
        subplot(1,3,TRIAL_FILTER)
        %get the distance and angle
        [~, DX, ~, OoBLen, flagOoB] = getDistAngle(dat, TRIAL_FILTER);
        Dists = getDistforOoB(DX, OoBLen, flagOoB);
        
        numSubjs = size(Dists,2);
        cm = jet(numSubjs); %pick numSubjs colors from 'jet'
        for subidx = 1:numSubjs
            dist = Dists{subidx};
            numtrials = size(dist,1);
            for nt = 1:numtrials
                pt = plot([1, 2], dist(nt,:), '-o', "Color",cm(subidx,:),'linewidth',2);
                %transparancy
                pt.Color = [pt.Color 0.3];
                hold on
            end
        
        end
        
        ylabel('Distance (meters)', 'FontSize', 15);
        xlim([0.5,2.5]);
        xticks([1 2]);
        xticklabels({'Actual','Reconstructed'});
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',10);
        title(trialname(TRIAL_FILTER), 'FontSize', 15)
    end
    

    sgtitle("GroupName: "+datname, 'FontSize', 20);
    
    exportgraphics(f,"C:/Users/Zilong/Desktop/path integration model/Andrea's matlab code/" + ...
    "PathIntegrationDataEstimationModelling/Output/DataVisualization/ReconstructDistDistribution/"+datname+".png",'Resolution',300);

    close(f);
end 


%% function to get distance 3
function Dists=getDistforOoB(DX, OoBLen, flagOoB)
%%plot angle 3 distibution
% DX is a cell structure containing the segment of each trial
% flagOoB is the flag of Out-of-Boundary trials, with 0 non-OoB trials 1 those OoB trials

subjectSize = size(DX,2);

Dists = cell(1,subjectSize);

for subj=1:subjectSize
    subjDX = DX{subj};
    subjOoBLen = OoBLen{subj};
    subjflagOoB = flagOoB{subj};

    sampleSize = size(subjDX,2);

    dist = [];

    for tr_id=1:sampleSize
        if subjflagOoB(tr_id)==0
            %non-OoB trial
            continue;
        else
            %OoB trial
            actual=subjOoBLen(tr_id); %actual length
            reconstructed=subjDX{tr_id}(3);  %reconstructed length
            dist=[dist;[actual,reconstructed]];
        end
    end
    Dists{subj}=dist;
end
end
