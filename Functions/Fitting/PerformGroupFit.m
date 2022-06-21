function Results = PerformGroupFit(GroupData, config)
%%Fit parameters for a single group
% GroupData is the extracted data containing the tracked positions of
% the participants as well as the their responses
% nsamples is the number of random samples to be taken when modelling the
% final position

TRIAL_FILTER = config.TrialFilter;          %load configurations necessary for the script
subjectNum = size(GroupData.FlagPos,2);     % Calculating sample size

%% Initialize empty cell for storing data
X               =       cell(1, subjectNum);% Actual positions
DX              =       cell(1, subjectNum);% Distances between subsequent points - segments li
THETADX         =       cell(1, subjectNum);% These are the angles between two subsequent segments. The angle indicate the rotation from the first segment towards the second, so it s the outer angle of the triangle.
segments        =       cell(1, subjectNum);
correctReDist   =       cell(1, subjectNum);
correctReAngle  =       cell(1, subjectNum);
DistErr         =       cell(1, subjectNum);
AngleErr        =       cell(1, subjectNum);
ProjSpeedL1     =       cell(1, subjectNum);% projected speed within detected start-to-end time window at leg 1
ProjSpeedL2     =       cell(1, subjectNum);% projected speed within detected start-to-end time window at leg 2
L1Dur           =       cell(1, subjectNum);% walking duration at leg 1
L2Dur           =       cell(1, subjectNum);% walking duration at leg2
StandingDur     =       cell(1, subjectNum);% standing duration at cone2
flagpos         =       cell(1, subjectNum);% flagpos
flagOoB         =       cell(1, subjectNum);% OoB flag
GroupParameters =       cell(1, subjectNum);% Output value
IC              =       cell(1, subjectNum);  

%% help function to calculate the angle between two vector
anglebetween = @(va,vb) atan2d(va(:,1).*vb(:,2) - va(:,2).*vb(:,1), va(:,1).*vb(:,1) + va(:,2).*vb(:,2));

%%
for j = 1:subjectNum

    %filter out participants who did short walking, happend in Young and HealthyOld, those out of distribution 
    if ismember(j, GroupData.BadPptIdxs)
        % set results to nan for later processing
        GroupParameters{j}  =   NaN(config.NumParams,1);
        IC{j}.aic           =   nan;
        IC{j}.bic           =   nan;
        IC{j}.negll         =   nan;
        IC{j}.likelihood    =   nan;
        flagOoB{j}          =   [];
        disp(['%%%%%%%%%%%%%%% Skipping PARTICIPANT ' num2str(j) ' ---- because they did a short walk%%%%%%%%%%%%%%%']);
        continue
    end
    
    if(TRIAL_FILTER == 0)
        %merge all three conditions
        flagpos{j}  = GroupData.FlagPos{j};
        BadExecutionTrials  = GroupData.Reconstructed{j}.BadExecution;
        realReturnAngles    = GroupData.Reconstructed{j}.RealReturnAngle;
        TrialNum    = size(GroupData.TrigPos{j},1);
        for idx = 1:TrialNum
            %If not out of bound or out of bound data is not present then take the trigpos
            if(GroupData.CondTable{j}.OutOfBound(idx) == 0 | isnan(GroupData.OutOfBoundPos{1,j}{idx}(1,1)))
                finalpos{j}(idx) = GroupData.TrigPos{j}(idx);
                flagOoB{j}(idx) = 0; %OoB flag is 0, i.e., not OoB trial
            elseif(GroupData.CondTable{j}.OutOfBound(idx) == 1)
                finalpos{j}(idx) = GroupData.ReconstructedOOB{j}.ReconstructedOoB(idx);
                flagOoB{j}(idx) = 1; %OoB flag is 0, i.e., it is OoB trial
            end
        end   
        Idx_GoodTrials      = BadExecutionTrials == 0;
        flagpos{j}          = flagpos{j}(Idx_GoodTrials);
        realReturnAngles    = realReturnAngles(Idx_GoodTrials);
        finalpos{j}         = finalpos{j}(Idx_GoodTrials);
        flagOoB{j}          = flagOoB{j}(Idx_GoodTrials); 
    
        leg1_duration       = GroupData.Reconstructed{j}.T_L1(Idx_GoodTrials);          %filter the duration of subject j at leg 1
        leg2_duration       = GroupData.Reconstructed{j}.T_L2(Idx_GoodTrials);          %filter the duration of subject j at leg 2
        standing_duration   = GroupData.Reconstructed{j}.T_Standing(Idx_GoodTrials);    %filter the duration of subject j standing at cone2 
        TrackedL1           = GroupData.TrackedL1{j}(Idx_GoodTrials);
        TrackedL2           = GroupData.TrackedL2{j}(Idx_GoodTrials);
    else
        %%filter trails based on the TRIAL_FILTER, with 1 no change, 2 no distal cues, 3 no optic flow
        Idx_Cond            = GroupData.CondTable{j}.Condition == TRIAL_FILTER; 
        flagpos{j}          = GroupData.FlagPos{j}(Idx_Cond);
        BadExecutionTrials  = GroupData.Reconstructed{j}.BadExecution(Idx_Cond);
        realReturnAngles    = GroupData.Reconstructed{j}.RealReturnAngle(Idx_Cond);
    
        tempCnt             = 1;
        TrialNum            = size(GroupData.TrigPos{j},1);
    
        %get the final position and OoB flag 
        for idx = 1:TrialNum %for each trial, if belongs to TRIAL_FILTER, go into
            if(GroupData.CondTable{j}.Condition(idx) == TRIAL_FILTER)
                %If not out of bound or out of bound data is not present then take the trigpos
                if(GroupData.CondTable{j}.OutOfBound(idx) == 0 | isnan(GroupData.OutOfBoundPos{1,j}{idx}(1,1)))
                    finalpos{j}(tempCnt) = GroupData.TrigPos{j}(idx);
                    flagOoB{j}(tempCnt) = 0; %OoB flag is 0, i.e., not OoB trial
                elseif(GroupData.CondTable{j}.OutOfBound(idx) == 1)
                    finalpos{j}(tempCnt) = GroupData.ReconstructedOOB{j}.ReconstructedOoB(idx);
                    flagOoB{j}(tempCnt) = 1; %OoB flag is 1, i.e., it is OoB trial
                end
                tempCnt = tempCnt + 1;
            end
        end
        
        Idx_GoodTrials      = BadExecutionTrials == 0;
        flagpos{j}          = flagpos{j}(Idx_GoodTrials);
        realReturnAngles    = realReturnAngles(Idx_GoodTrials);
        finalpos{j}         = finalpos{j}(Idx_GoodTrials);
        flagOoB{j}          = flagOoB{j}(Idx_GoodTrials); 
    
        leg1_duration       = GroupData.Reconstructed{j}.T_L1(Idx_Cond);          %filter the duration of subject j at leg 1
        leg1_duration       = leg1_duration(Idx_GoodTrials);  
    
        leg2_duration       = GroupData.Reconstructed{j}.T_L2(Idx_Cond);          %filter the duration of subject j at leg 2
        leg2_duration       = leg2_duration(Idx_GoodTrials);
    
        standing_duration   = GroupData.Reconstructed{j}.T_Standing(Idx_Cond);    %filter the duration of subject j standing at cone2 
        standing_duration   = standing_duration(Idx_GoodTrials);
    
        TrackedL1           = GroupData.TrackedL1{j}(Idx_Cond);
        TrackedL1           = TrackedL1(Idx_GoodTrials);
    
        TrackedL2           = GroupData.TrackedL2{j}(Idx_Cond);
        TrackedL2           = TrackedL2(Idx_GoodTrials);
    end
    
    %get X, segments, DX etc....
    for tr = 1:length(flagpos{j})
        
        X{j}{tr}         =      [flagpos{j}{tr}(:,[1,3]);finalpos{j}{tr}([1,3])];
        segments{j}{tr}  =      X{j}{tr}(2:end,:) - X{j}{tr}(1:end-1,:);
        DX{j}{tr}        =      sqrt(sum(segments{j}{tr}.^2,2));        
     
        outer_rad        =      [0; anglebetween(segments{j}{tr}(1:end-1,:), segments{j}{tr}(2:end,:))];
        outer_rad(3,1)   =      realReturnAngles(tr);     
        THETADX{j}{tr}   =      deg2rad(outer_rad);%wrap the angle into (0,2pi)
        
%         %extract correct return distance and angle
%         if flagOoB{j}(tr) == 0
%             x_2 = X{j}{tr}(3,:);
%             correctReDist{j}{tr} = sqrt(sum(x_2.^2));
% 
%             p1 = X{j}{tr}(1,:);
%             p2 = X{j}{tr}(2,:);
%             p3 = X{j}{tr}(3,:);
%             vec1 = p3-p2; 
%             vec2 = p1-p3;
%             correctReAngle{j}{tr} = anglebetween(vec1, vec2); 
% 
%             %calculate distance error and angular error
%             DistErr{j}{tr} = correctReDist{j}{tr}-DX{j}{tr}(3);
%             AngleErr{j}{tr} = correctReAngle{j}{tr}-realReturnAngles(tr);
%         else 
%             DistErr{j}{tr} = nan;
%             AngleErr{j}{tr} = nan;
%         end

        %extract correct return distance and angle
        if flagOoB{j}(tr) == 0
            x_2 = X{j}{tr}(3,:);
            correctReDist{j}{tr} = sqrt(sum(x_2.^2));
            DistErr{j}{tr} = correctReDist{j}{tr}-DX{j}{tr}(3);
        else 
            DistErr{j}{tr} = nan;
        end

        p1 = X{j}{tr}(1,:);
        p2 = X{j}{tr}(2,:);
        p3 = X{j}{tr}(3,:);
        vec1 = p3-p2; 
        vec2 = p1-p3;
        correctReAngle{j}{tr} = anglebetween(vec1, vec2); 
        AngleErr{j}{tr} = correctReAngle{j}{tr}-realReturnAngles(tr);

        %extract the projected speed information along with the time information on outbound path
        L1_Vel_proj             =       TrackedL1{tr}.Vel_proj;
        L1_Time                 =       TrackedL1{tr}.Time;
        L1_Filtered_Vel_proj    =       TrackedL1{tr}.Filtered_Vel_proj;
        L1_Vel_proj_selected    =       L1_Vel_proj(L1_Filtered_Vel_proj);
        L1_Time_selected        =       L1_Time(L1_Filtered_Vel_proj);
        ProjSpeedL1{j}{1,tr}    =       L1_Time_selected;  %selected time in a start-to-end range
        ProjSpeedL1{j}{2,tr}    =       L1_Vel_proj_selected;  %selected speed in a start-to-end range

        L2_Vel_proj             =       TrackedL2{tr}.Vel_proj;
        L2_Time                 =       TrackedL2{tr}.Time;
        L2_Filtered_Vel_proj    =       TrackedL2{tr}.Filtered_Vel_proj;
        L2_Vel_proj_selected    =       L2_Vel_proj(L2_Filtered_Vel_proj);
        L2_Time_selected        =       L2_Time(L2_Filtered_Vel_proj);
        ProjSpeedL2{j}{1, tr}   =       L2_Time_selected;  %selected time in a start-to-end range
        ProjSpeedL2{j}{2, tr}   =       L2_Vel_proj_selected; %selected speed in a start-to-end range
        
        L1Dur{j}{tr}            =       leg1_duration(tr);      %extract the walking duration at leg 1 
        L2Dur{j}{tr}            =       leg2_duration(tr);      %extract the walking duration at leg 2
        StandingDur{j}{tr}      =       standing_duration(tr);  %extract the standing duration at cone2
    end


    %% put all of the things we need into a struct for sending to FitData
    Input.DX                 =   DX{j};
    Input.THETADX            =   THETADX{j};
    Input.X                  =   X{j};
    Input.flagOoB            =   flagOoB{j};
    Input.ProjSpeedL1        =   ProjSpeedL1{j};
    Input.ProjSpeedL2        =   ProjSpeedL2{j};
    Input.L1Dur              =   L1Dur{j};
    Input.L2Dur              =   L2Dur{j};
    Input.StandingDur        =   StandingDur{j};

    if length(flagpos{j}) < config.NumParams
        %lack of trials, skip estimation
        disp("%%%%%%%%%%%%%%% Skipping participant " + num2str(j) + ...
            ", because only "+ length(flagpos{j}) + ...
            " datapoints available for parameter estimation%%%%%%%%%%%%%%%\n");
        % set results to nan for later processing
        GroupParameters{j}  =   NaN(config.NumParams,1);
        IC{j}.aic           =   nan;
        IC{j}.bic           =   nan;
        IC{j}.negll         =   nan;
        IC{j}.likelihood    =   nan;
        %flagOoB{j}          =   [];
        continue;
    end

    %% Do the data fitting
    disp(['%%%%%%%%%%%%%%% STARTING FIT PER PARTICIPANT ' num2str(j) ' %%%%%%%%%%%%%%%']);
    [GroupParameters{j}, IC{j}] = FitData(Input, config);
end

%%
%Transforming the fitted parameters from cell to array
[~, rows]       = size(GroupParameters);
[cols,~]        = size(GroupParameters{1});
estimatedParams = zeros(rows,cols);
for i=1:rows
   estimatedParams(i,:) = GroupParameters{i}';
end

%put all results into a matlab structure
Results.estimatedParams =   estimatedParams;
Results.X               =   X;
Results.DX              =   DX;
Results.THETADX         =   THETADX;
Results.IC              =   IC;
Results.flagOoB         =   flagOoB;
Results.DistErr         =   DistErr;
Results.AngleErr        =   AngleErr;

end