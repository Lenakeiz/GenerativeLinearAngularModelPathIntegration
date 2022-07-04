%% Preparing the data
VAM

%% Setting colors for using in plots
ColorPattern;
resultfolder = savefolder + "PaperFigs/Fig1A";
config.ResultFolder = resultfolder;

%%
[MCIPosNormedDist, MCIPosNormedAngle] = getNormalizedReturnDistanceandAngle(MCIPos.X, MCIPos.DX, MCIPos.Theta)

%% get the normalized return distance and angles for all particpants
% each value is a mean value for per participant
function [NormedDist, NormedAngle] = getNormalizedReturnDistanceandAngle(AllX, AllDX, AllTheta)
    numConds = length(AllX);
    NormedDist = [];
    NormedAngle = [];

    for TRIAL_FILTER=1:numConds
        X = AllX{TRIAL_FILTER};
        DX = AllDX{TRIAL_FILTER};
        Theta = AllTheta{TRIAL_FILTER};

        subjectSize = size(X,2);
        NormedDistAllSubjs = zeros(1,subjectSize);
        NormedAngleAllSubjs = zeros(1,subjectSize);

        for subj=1:subjectSize
            subjX = X{subj};
            subjDX = DX{subj};
            subjTheta = Theta{subj};

            sampleSize = size(subjX,2);
            NormedDistPerSubj = zeros(1,sampleSize);
            NormedAnglePerSubj = zeros(1,sampleSize);
            for tr_id=1:sampleSize
                %normalized return distance
                %actual return distance
                l3 = subjDX{tr_id}(3);
                %correct return distance
                p3 = subjX{tr_id}(3,:);
                correct_l3 = norm(p3);

                normalized_dist = l3/correct_l3;
                NormedDistPerSubj(tr_id) = normalized_dist;

                %normalized return angle
                %actual return angle
                theta3 = subjTheta{tr_id}(3);
                %correct return angle 
                p1 = subjX{tr_id}(1,:);
                p2 = subjX{tr_id}(2,:);
                p3 = subjX{tr_id}(3,:);
                vec1 = p3-p2; vec2 = p1-p3;
                correct_theta3=atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
                correct_theta3 = deg2rad(correct_theta3);

                normalized_angle = theta3/correct_theta3;
                NormedAnglePerSubj(tr_id) = normalized_angle;
            end
            NormedDistAllSubjs(subj)=mean(NormedDistPerSubj);
            NormedAngleAllSubjs(subj)=mean(NormedAnglePerSubj);
        end
        NormedDist = [NormedDist;NormedDistAllSubjs]; %dim=3*NumSubjects
        NormedAngle = [NormedAngle;NormedAngleAllSubjs]; %dim=3*NumSubjects
    end
end
