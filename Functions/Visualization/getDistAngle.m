function [X, DX, THETADX, OoBLen, flagOoB] = getDistAngle(GroupData, TRIAL_FILTER)
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
% OoB flag
flagOoB = cell(1, sampleSize);
% OoB length
OoBLen = cell(1, sampleSize);

anglebetween = @(va,vb) atan2d(va(:,1).*vb(:,2) - va(:,2).*vb(:,1), va(:,1).*vb(:,1) + va(:,2).*vb(:,2));

for j = 1:sampleSize

    if(TRIAL_FILTER == 0)
        %processing the data from all conditions
        flagpos{j}  = GroupData.FlagPos{j};
        OoBLen{j} = GroupData.Errors{j}.OoBLength;
        for idx = 1:size(GroupData.TrigPos{j},1)
            %If not out of bound or out of bound data is not present then take the trigpos
            if(GroupData.CondTable{j}.OutOfBound(idx) == 0 | isnan(GroupData.OutOfBoundPos{1,j}{idx}(1,1)))
                finalpos{j,1}(idx,1) = GroupData.TrigPos{j}(idx);
                flagOoB{j}(idx) = 0; %OoB flag is 0, i.e., not OoB trial
            elseif(GroupData.CondTable{j}.OutOfBound(idx) == 1)
                finalpos{j,1}(idx,1) = GroupData.ReconstructedOOB{j}.ReconstructedOoB(idx);
                flagOoB{j}(idx) = 1; %OoB flag is 0, i.e., it is OoB trial
            end
        end
    else
        %processing data according to "TRIAL_FILTER", with 1 no change, 2 no distal cues, 3 no optic flow
        flagpos{j}  = GroupData.FlagPos{j}(GroupData.CondTable{1,j}.Condition == TRIAL_FILTER);
        OoBLen{j} = GroupData.Errors{j}.OoBLength(GroupData.CondTable{1,j}.Condition == TRIAL_FILTER);
        tempCnt = 1;
        for idx = 1:size(GroupData.TrigPos{j},1)
            if(GroupData.CondTable{1,j}.Condition(idx) == TRIAL_FILTER)
                %If not out of bound or out of bound data is not present then take the trigpos
                if(GroupData.CondTable{j}.OutOfBound(idx) == 0 | isnan(GroupData.OutOfBoundPos{1,j}{idx}(1,1)))
                    finalpos{j,1}(tempCnt,1) = GroupData.TrigPos{j}(idx);
                    flagOoB{j}(tempCnt) = 0; %OoB flag is 0, i.e., not OoB trial
                elseif(GroupData.CondTable{j}.OutOfBound(idx) == 1)
                    finalpos{j,1}(tempCnt,1) = GroupData.ReconstructedOOB{j}.ReconstructedOoB(idx);
                    flagOoB{j}(tempCnt) = 1; %OoB flag is 0, i.e., it is OoB trial
                end
                tempCnt = tempCnt + 1;
            end
        end
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