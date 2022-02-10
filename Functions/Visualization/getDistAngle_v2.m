function [X, DX, THETADX] = getDistAngle_v2(GroupData, TRIAL_FILTER)
%% extract detailed distance and angle information of the group data
% GroupData is the extracted data containing the tracked positions of
% the participants as well as the their responses
% nsamples is the number of random samples to be taken when modelling the
% final position

% Calculating sample size
sampleSize = size(GroupData.FlagPos,2);

% Actual positions
X = cell(1, sampleSize);
% Distances between subsequent points - segments li
DX = cell(1, sampleSize);
% These are the angles between two subsequent segments. The angle indicate
% the rotation from the first segment towards the second, so it s the outer
% angle of the triangle.
THETADX = cell(1, sampleSize);

anglebetween = @(va,vb) atan2d(va(:,1).*vb(:,2) - va(:,2).*vb(:,1), va(:,1).*vb(:,1) + va(:,2).*vb(:,2));

for j = 1:sampleSize

    if(TRIAL_FILTER == 0)
        flagpos{j}  = GroupData.FlagPos{j}(GroupData.Cond{j}(:,3) == 0);
        finalpos{j} = GroupData.TrigPos{j}(GroupData.Cond{j}(:,3) == 0);
    else
        flagpos{j}  = GroupData.FlagPos{j}(GroupData.Cond{j}(:,3) == 0 & GroupData.CondTable{1,j}.Condition == TRIAL_FILTER);
        finalpos{j} = GroupData.TrigPos{j}(GroupData.Cond{j}(:,3) == 0 & GroupData.CondTable{1,j}.Condition == TRIAL_FILTER);
    end


    % Create structures to save the data for each trial
    DX{j}       = cell(1,length(flagpos{j}));
    THETADX{j}  = cell(1,length(flagpos{j}));
    X{j}        = cell(1,length(flagpos{j}));
    
    % Temporary
    segments = cell(1, length(sampleSize));

    for tr = 1:length(flagpos{j})
        
        X{j}{tr}         = [flagpos{j}{tr}(:,[1,3]);finalpos{j}{tr}([1,3])];
        %Each segment is just the difference between the trials of the 
        segments{j}{tr}  = X{j}{tr}(2:end,:) - X{j}{tr}(1:end-1,:);
        DX{j}{tr}        = sqrt(sum(segments{j}{tr}.^2,2));        
        
        %THETADX{j}{tr} = ([0; anglebetween(segments{j}{tr}(1:end-1,:), segments{j}{tr}(2:end,:))]);
        outer_rad = deg2rad([0; anglebetween(segments{j}{tr}(1:end-1,:), segments{j}{tr}(2:end,:))]);

        %MODIFICATION!ATTENTION!:wrap the angle to [0-2pi)
        %THETADX{j}{tr} = outer_rad;
        THETADX{j}{tr} = mod(outer_rad, 2*pi);
        
    end

end
end