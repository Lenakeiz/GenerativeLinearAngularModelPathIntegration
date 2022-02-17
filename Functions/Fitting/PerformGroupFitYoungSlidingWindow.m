function Results = PerformGroupFitYoungSlidingWindow(GroupData, windowsize, config)
%% Fit Theta parameters for a single group
% GroupData is the extracted data containing the tracked positions of
% the participants as well as the their responses
% windowsize, =0.7 means using the 70% data 

%load configurations necessary for the script
TRIAL_FILTER = config.TrialFilter;
numParams = config.NumParams;

%obtain all the square sum of distances, l_1^2+l_2^2
SSOD = getAllSquareSumOfDistance(GroupData,TRIAL_FILTER);

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

SelectedOoBLen = cell(1, sampleSize);
SelectedflagOoB = cell(1, sampleSize);

% Output value
GroupParameters = cell(1, sampleSize);
BIC = cell(1, sampleSize);

anglebetween = @(va,vb) atan2d(va(:,1).*vb(:,2) - va(:,2).*vb(:,1), va(:,1).*vb(:,1) + va(:,2).*vb(:,2));

numdat = length(SSOD);

interval = 0.05;
Results = cell(1, floor((1-windowsize)/interval)+1);
index = 1;
for lowpercent=0:interval:1-windowsize
    uppercent = lowpercent+windowsize;
    %find lower bound and upper bound of the sliding window
    sortSSOD = sort(SSOD, 'ascend');
    lowidx = ceil(numdat*lowpercent+eps); %add eps to make sure we won't get index 0 
    lowbound= sortSSOD(lowidx);
    upidx = floor(numdat*uppercent); 
    upbound = sortSSOD(upidx);

    %do the fitting
    for j = 1:sampleSize
        
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
    
        if length(flagpos{j}) < numParams
            disp("%%%%%%%%%%%%%%% WARNING: " + length(flagpos{j}) + " datapoints for the current participant %%%%%%%%%%%%%%%\n");
        end
    
        if isempty(flagpos{j})
            disp("%%%%%%%%%%%%%%% Skipping participant " + num2str(j) + " because no datapoints available %%%%%%%%%%%%%%%");
            % NUmber of parameters in return matrix is fixed
            GroupParameters{j} = NaN(8,1);
            BIC{j} = nan;
            flagOoB{j} = [];
            continue;
        end
        
        % Temporary
        segments = cell(1, length(sampleSize));
        
        temptr = 1;
        for tr = 1:length(flagpos{j})
            
            tempX = flagpos{j}{tr}(:,[1,3]);
            tempsegments = tempX(2:end,:) - tempX(1:end-1,:);
            tempssod = sum(tempsegments.^2,[1,2]);

            %if the current outbound square sum  is in the selected range
            if tempssod>=lowbound && tempssod<=upbound
                X{j}{temptr}         = [flagpos{j}{tr}(:,[1,3]);finalpos{j}{tr}([1,3])];
                %Each segment is just the difference between the trials of the 
                segments{j}{temptr}  = X{j}{temptr}(2:end,:) - X{j}{temptr}(1:end-1,:);
                DX{j}{temptr}        = sqrt(sum(segments{j}{temptr}.^2,2)); 

                SelectedOoBLen{j}(temptr) = OoBLen{j}(tr);
                SelectedflagOoB{j}(temptr) = flagOoB{j}(tr);
                
                outer_rad = deg2rad([0; anglebetween(segments{j}{temptr}(1:end-1,:), segments{j}{temptr}(2:end,:))]);
        
                %MODIFICATION!ATTENTION!:wrap the angle to [0-2pi)
                THETADX{j}{temptr} = mod(outer_rad, 2*pi);

                temptr = temptr+1;
            else
                continue;
            end     
        end
        
        if isempty(X{j})
            %no enough trials, set the parameter fitting of the participant
            %to Nan
            GroupParameters{j} = nan(8,1);
        else
            disp(['%%%%%%%%%%%%%%% STARTING FIT PER PARTICIPANT ' num2str(j) ' %%%%%%%%%%%%%%%']);
            % Do the data fitting
            [GroupParameters{j}, BIC{j}] = FitData(DX{j},THETADX{j},X{j},SelectedOoBLen{j},SelectedflagOoB{j},config);
        end
    
    end

    %Transforming the fitted parameters to array
    [~, rows] = size(GroupParameters);
    [cols,~] = size(GroupParameters{1});
    group_results = zeros(rows,cols);
    for i=1:rows
       group_results(i,:) = GroupParameters{i}';
    end

    Results{index} = group_results;
    index = index+1;

end 

end

function SSOD = getAllSquareSumOfDistance(GroupData,TRIAL_FILTER) 
    % Calculating sample size
    subjSize = size(GroupData.FlagPos,2);
    SSOD = [];
    for subj = 1:subjSize
        flagpos{subj}  = GroupData.FlagPos{subj}(GroupData.CondTable{1,subj}.Condition == TRIAL_FILTER);
        for trialidx = 1:length(flagpos{subj})
            X = flagpos{subj}{trialidx}(:,[1,3]);
            segment  = X(2:end,:) - X(1:end-1,:);
            ssod = sum(segment.^2,[1,2]);
            SSOD = [SSOD,ssod];
        end
    end

end