function Results = PerformGroupFit(GroupData, config)
%%Fit parameters for a single group
% GroupData is the extracted data containing the tracked positions of
% the participants as well as the their responses
% nsamples is the number of random samples to be taken when modelling the
% final position

%load configurations necessary for the script
TRIAL_FILTER = config.TrialFilter;
numParams = config.NumParams;
Model_Name = config.ModelName;

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
% projected speed within detected start-to-end time window
ProjSpeedL1 = cell(1, sampleSize);
ProjSpeedL2 = cell(1, sampleSize);

% flagpos
flagpos = cell(1, sampleSize);
% OoB flag
flagOoB = cell(1, sampleSize);
% OoB length
OoBLen = cell(1, sampleSize);

% Output value
GroupParameters = cell(1, sampleSize);
IC = cell(1, sampleSize);  

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

    if length(flagpos{j}) < numParams
        disp("%%%%%%%%%%%%%%% WARNING: " + length(flagpos{j}) + " datapoints for the current participant %%%%%%%%%%%%%%%\n");
    end

    if isempty(flagpos{j}) %this can be empty if we remove the OoB trials. Empty coz all the trials are OoB trials
        disp("%%%%%%%%%%%%%%% Skipping participant " + num2str(j) + " because no datapoints available %%%%%%%%%%%%%%%");
        % NUmber of parameters in return matrix is fixed
        if Model_Name == "G1G2"%G1G2 model has 8 parameters in total
            GroupParameters{j} = NaN(8,1);
        else %gamma model has 7 parameters in total
            GroupParameters{j} = NaN(7,1);
        end
        IC{j} = nan;
        flagOoB{j} = [];
        continue;
    end

    % Create structures to save the data for each trial
    DX{j}       = cell(1,length(flagpos{j}));
    THETADX{j}  = cell(1,length(flagpos{j}));
    X{j}        = cell(1,length(flagpos{j}));
    ProjSpeedL1{j}= cell(2,length(flagpos{j}));
    ProjSpeedL2{j}= cell(2,length(flagpos{j}));

    % Temporary
    segments = cell(1, length(sampleSize));

    for tr = 1:length(flagpos{j})
        
        X{j}{tr}         = [flagpos{j}{tr}(:,[1,3]);finalpos{j}{tr}([1,3])];
        %Each segment is just the difference between the trials of the 
        segments{j}{tr}  = X{j}{tr}(2:end,:) - X{j}{tr}(1:end-1,:);
        DX{j}{tr}        = sqrt(sum(segments{j}{tr}.^2,2));        
        
        %THETADX{j}{tr} = ([0; anglebetween(segments{j}{tr}(1:end-1,:), segments{j}{tr}(2:end,:))]);
        outer_rad = deg2rad([0; anglebetween(segments{j}{tr}(1:end-1,:), segments{j}{tr}(2:end,:))]);

        %wrap the angle into (0,2pi)
        THETADX{j}{tr} = mod(outer_rad, 2*pi);

        %
        L1_Vel_proj = GroupData.TrackedL1{j}{tr}.Vel_proj;
        L1_Time = GroupData.TrackedL1{j}{tr}.Time;
        L1_Filtered_Vel_proj = GroupData.TrackedL1{j}{tr}.Filtered_Vel_proj;
        L1_Vel_proj_selected = L1_Vel_proj(L1_Filtered_Vel_proj);
        L1_Time_selected = L1_Time(L1_Filtered_Vel_proj);
        ProjSpeedL1{j}{1,tr} = L1_Time_selected;  %selected time in a start-to-end range
        ProjSpeedL1{j}{2,tr} = L1_Vel_proj_selected;  %selected speed in a start-to-end range

        L2_Vel_proj = GroupData.TrackedL2{j}{tr}.Vel_proj;
        L2_Time = GroupData.TrackedL2{j}{tr}.Time;
        L2_Filtered_Vel_proj = GroupData.TrackedL2{j}{tr}.Filtered_Vel_proj;
        L2_Vel_proj_selected = L2_Vel_proj(L2_Filtered_Vel_proj);
        L2_Time_selected = L2_Time(L2_Filtered_Vel_proj);
        ProjSpeedL2{j}{1, tr} = L2_Time_selected;  %selected time in a start-to-end range
        ProjSpeedL2{j}{2, tr} = L2_Vel_proj_selected; %selected speed in a start-to-end range
        
    end

    % Do the data fitting
    disp(['%%%%%%%%%%%%%%% STARTING FIT PER PARTICIPANT ' num2str(j) ' %%%%%%%%%%%%%%%']);
    [GroupParameters{j}, IC{j}] = FitData(DX{j}, THETADX{j}, X{j}, ProjSpeedL1{j}, ProjSpeedL2{j}, config);
end

%Transforming the fitted parameters to array
[~, rows] = size(GroupParameters);
[cols,~] = size(GroupParameters{1});
estimatedParams = zeros(rows,cols);
for i=1:rows
   estimatedParams(i,:) = GroupParameters{i}';
end

%put all results into a matlab structure
Results.estimatedParams = estimatedParams;
Results.X = X;
Results.DX = DX;
Results.THETADX = THETADX;
Results.IC = IC;
Results.flagOoB = flagOoB;

end