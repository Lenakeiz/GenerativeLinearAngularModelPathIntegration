function [outGroup] = CalculateTrackingPath(Group, config)
%% CalculateTrackingPath
% Andrea Castegnaro, UCL, 2022 uceeaca@ucl.ac.uk
% Calculate trial information based on tracking data.
% Despite having a separate output in our dataset (i.e. cones locations,
% cones triggered position, cones spawning times, triggered position, out
% of bound location on incorrect trials) we wanted to analyse the tracking
% data saved for each trial. We reconstructed from the tracking data the
% following quantities. Speed along segments l1 and l2, body angular
% rotation at cone 3, reconstructed walked lengths of l1 and l2, time spent
% on l1 and l2, time spent standing at cone 2.
% We used the body angular rotation at cone 3 to detect trials that needed
% a visual check, to see if participants completed the triangle using the
% short or long rotation, to account for the correct angular execution
% error.
% ===================================================================================

pSize = size(Group.Path,2);

Group.TrackedL1 = {};
Group.TrackedL2 = {};
Group.Reconstructed = {};
Group.BadPptIdxs = [];
Group.overRotations = [];

% Anonymous function to calculate angle between two vectors.
anglebetween = @(v,w) atan2d(w(:,2).*v(:,1) - v(:,2).*w(:,1), v(:,1).*w(:,1) + v(:,2).*w(:,2));

% Delta T set to integrate the angular information from the tracking data
dtAngular = config.TrackedInboundAngularDeltaT;

for pId = 1:pSize

    trialSize = size(Group.Path{1,pId},1);

    Group.TrackedL1{1,pId} = {};
    Group.TrackedL2{1,pId} = {};

    % Preparing output structures
    InboundAngularRotation = [];
    SignInboundBodyRotation = [];
    InferredReturnAngle      = [];
    SignInferredReturnAngle = [];
    OobReturnAngle = [];
    SignOobReturnAngle = [];
    RealReturnAngle = [];

    participantRecontructedQuantities = [];

    % Calculate the mean time between making cone 1 disappear and spawning cone
    % two (same between 2 and 3) animations
    % It should be approximately 2 seconds
    DiffTrigSpawn12 = cell2mat(Group.FlagSpawnTime{1,pId}(:,2)) - cell2mat(Group.FlagTrigTimes{1,pId}(:,1));
    DiffTrigSpawn23 = cell2mat(Group.FlagSpawnTime{1,pId}(:,3)) - cell2mat(Group.FlagTrigTimes{1,pId}(:,2));

    for trialId = 1:trialSize

        Cone_pos    = Group.FlagPos{1,pId}{trialId,1};
        Trig_pos    = Group.TrigPos{1,pId}{trialId,1};
        Tracked_pos = Group.Path{1,pId}{trialId,1};
        Tracked_pos = array2table(Tracked_pos,"VariableNames",{'Time' 'Pos_X' 'Pos_Y' 'Pos_Z' 'Forward_X' 'Forward_Y' 'Forward_Z'});

        % Getting Out of Bound (OoB) location in OoB trials 
        oob_return_angle = nan;
        outofbound = Group.CondTable{1,pId}.OutOfBound(trialId);
        if (outofbound == 1)
            OoB_pos = Group.OutOfBoundPos{1,pId}{trialId,1};
            
            % If CalculateTrackingPath is run after TransformPaths then we 
            % can use the reconstructed OoB location. See details in
            % Reconstructed OoB location have a distance from cone 3 set as
            % the average distance calculated for all good trials in that
            % condition. This is useful for the behavioural analyisis where
            % distance will not weight on the mean calculation. In the
            % model fitting the distance is not considered in OoB trials,
            % but only the angular rotation.
            if(isfield(Group,'ReconstructedOOB'))
                OoB_pos = Group.ReconstructedOOB{1,pId}.ReconstructedOoB{trialId,1};
            end
        end

        %% Calculating angular rotations from tracking data
        % Filtering tracked position after reaching cone 3
        Angular_pos = Tracked_pos(Tracked_pos.Time > Group.FlagTrigTimes{1,pId}{trialId,3} - 2,:);
        % Getting only Forward X - Z components
        Angular_pos = removevars(Angular_pos,{'Pos_X' 'Pos_Y' 'Pos_Z' 'Forward_Y'});
        % Getting a reading every dtAngular
        Angular_pos = Angular_pos(1:dtAngular:end,:);

        % Calculating the signed angle between vectors indexed by cone 2-3
        % and indexed by triggered position - cone 3 
        dir_23   = [(Cone_pos(3,1) - Cone_pos(2,1)) (Cone_pos(3,3) - Cone_pos(2,3))];
        dir_trig = [(Trig_pos(1,1) - Cone_pos(3,1)) (Trig_pos(1,3) - Cone_pos(3,3))];
        f_return_angle = anglebetween(dir_23,dir_trig);
        
        % Similar but if OoB trial we use the OoB location
        if (outofbound == 1)
            dir_oob = [(OoB_pos(1,1) - Cone_pos(3,1)) (OoB_pos(1,3) - Cone_pos(3,3))];
            oob_return_angle = anglebetween(dir_23,dir_oob);
        end

        % Calculating the body rotation for this trial by summing up
        % all of the signed angular dts
        tracking_size = height(Angular_pos);
        ang_rotation = 0;
        for trackId = 2:tracking_size
            prevDirection = [Angular_pos.Forward_X(trackId-1) Angular_pos.Forward_Z(trackId-1)];
            prevDirection = prevDirection/norm(prevDirection);
            currDirection = [Angular_pos.Forward_X(trackId) Angular_pos.Forward_Z(trackId)];
            currDirection = currDirection/norm(currDirection);
            ang_rotation = ang_rotation + anglebetween(prevDirection,currDirection);
        end

        % If the cumulative body rotation is higher (lower) than 180 (-180)
        % this means the participant has looked around more on the spot. We
        % flag this as an over rotation, because it will affect the
        % executing return angle by the participant. Note this could not be
        % detected from the simple response (OoB location) of the
        % participant in relation to cone 3 
        if(ang_rotation > 180 || ang_rotation < -180)
            disp(['Over rotation participant ', num2str(pId), ' trial ', num2str(trialId)]);
            Group.overRotations = [Group.overRotations;[pId trialId]];
        end

        % Based on tracked angular rotation we can calculate whether there
        % is a mismatch between the real angular body rotation and the
        % angular rotation inferred by triggered position. Let s calculate
        % the angular rotation sign
        sign_ang_rotation = sign(ang_rotation); % as calculated from the tracking data 
        sign_return_angle = sign(f_return_angle); % as calculated from the triggered location
        sign_oob_return_angle = sign(oob_return_angle); % as calculated from the OoB lcoation
 
        % Let s create temporary variables for the participant performed
        % angular rotation at cone 3. If OoB trial instead of triggered
        % location we use the inferred angular. We want to calculate the
        % angular rotation from the participant response (OoB location),
        % but we need to account for the fact participant rotated over the
        % short angle or long angle. We now know this information from the
        % calculation made above, so we are going to offset the rotation
        % based on this information.
        if(outofbound == 1)
            temp_ret_angle = oob_return_angle;
            tem_sign_ret_angle = sign_oob_return_angle;
        else
            temp_ret_angle = f_return_angle;
            tem_sign_ret_angle = sign_return_angle;
        end

        % Offseting the real angle based on our calculation on
        % partiticipant rotation at cone 3 from tracking data
        real_return_angle = temp_ret_angle;
        if(sign_ang_rotation * tem_sign_ret_angle == -1)
            % We found a mismatch between inferred rotation and actual
            % body rotation. In this case we are going to readjust the
            % angles. Positive is an anticlockwise rotation, while
            % negative is clockwise rotation
            if (tem_sign_ret_angle > 0)
                real_return_angle = real_return_angle - 360;
            elseif (tem_sign_ret_angle < 0)
                real_return_angle = real_return_angle + 360;
            end
        end

        InboundAngularRotation  = [InboundAngularRotation;ang_rotation]; % body rotation
        SignInboundBodyRotation = [SignInboundBodyRotation;sign_ang_rotation]; % body rotation sign 
        InferredReturnAngle     = [InferredReturnAngle;f_return_angle]; % Inferred from triggerered location
        SignInferredReturnAngle = [SignInferredReturnAngle;sign_return_angle]; % Inferred from triggerered pos sign
        OobReturnAngle          = [OobReturnAngle;oob_return_angle]; % Inferred from OoB location
        SignOobReturnAngle      = [SignOobReturnAngle;sign_oob_return_angle]; % Inferred from OoB location sign

        % From trig location or oob location depending on cases. This will
        % be used in our model.
        RealReturnAngle         = [RealReturnAngle;real_return_angle]; 

        %% Calculating speed and segment lengths (l1 l2) from tracking data
        % Note that segment lengths have been used only for sanity check
        % for the calculated speed and not input directly into the model
        RecontructedQuantities = nan(1,13);

        % Filtering tracking data by excluding locations before reaching
        % cone 1 
        Tracked_pos    = Tracked_pos(Tracked_pos.Time > Group.FlagTrigTimes{1,pId}{trialId,1},:);

        % Extracting tracking data for l1 segment (between cone 2 spawn
        % time and cone 2 reaching timestamps)
        Tracked_pos_L1 = Tracked_pos(Tracked_pos.Time > Group.FlagSpawnTime{1,pId}{trialId,2} &...
            Tracked_pos.Time < Group.FlagTrigTimes{1,pId}{trialId,2} + mean(DiffTrigSpawn12) + config.Speed.timeOffsetAfterFlagReach, :);

        % Extracting tracking data for l2 segment (between cone 3 spawn and
        % time and cone 3 reaching timestamps)
        Tracked_pos_L2 = Tracked_pos(Tracked_pos.Time > Group.FlagSpawnTime{1,pId}{trialId,3} &...
            Tracked_pos.Time < Group.FlagTrigTimes{1,pId}{trialId,3} + mean(DiffTrigSpawn23) + config.Speed.timeOffsetAfterFlagReach, :);

        % Calculating ideal direction of walking along the two segments
        s_dir_L1 = (Cone_pos(2,[1 3]) - Cone_pos(1,[1 3])) / norm(Cone_pos(2,[1 3]) - Cone_pos(1,[1 3]));
        s_dir_L2 = (Cone_pos(3,[1 3]) - Cone_pos(2,[1 3])) / norm(Cone_pos(3,[1 3]) - Cone_pos(2,[1 3]));

        %% Calculating a moving average velocity while walking on l1 segment
        l1Size = height(Tracked_pos_L1);
        Tracked_pos_L1.Vel_X = zeros(l1Size,1);
        Tracked_pos_L1.Vel_Z = zeros(l1Size,1);

        for iL = 2:l1Size

            Tracked_pos_L1.Vel_X(iL) = config.Speed.alpha * Tracked_pos_L1.Vel_X(iL-1) +...
                (1 - config.Speed.alpha) * (Tracked_pos_L1.Pos_X(iL) - Tracked_pos_L1.Pos_X(iL-1)) / (Tracked_pos_L1.Time(iL) - Tracked_pos_L1.Time(iL-1));

            Tracked_pos_L1.Vel_Z(iL) = config.Speed.alpha * Tracked_pos_L1.Vel_Z(iL-1) +...
                (1 - config.Speed.alpha) * (Tracked_pos_L1.Pos_Z(iL) - Tracked_pos_L1.Pos_Z(iL-1)) / (Tracked_pos_L1.Time(iL) - Tracked_pos_L1.Time(iL-1));

        end

        % Calculating projection of calculated velocity over ideal walking
        % direction
        Tracked_pos_L1.Vel_proj = dot([Tracked_pos_L1.Vel_X Tracked_pos_L1.Vel_Z],repmat(s_dir_L1,l1Size,1),2);
        % Gaussian filter to smooth data
        Tracked_pos_L1.Smoothed_Vel_proj = smoothdata(Tracked_pos_L1.Vel_proj,'gaussian',config.Speed.smoothWindow);

        % Using a speed threshold cutoff we can isolate from our
        % reconstruction only the data where the participant actually moved 
        filteredVel_proj = Tracked_pos_L1.Smoothed_Vel_proj > config.Speed.velocityCutoff;
        [~,currIdxsL1] = findWalkingSegment(filteredVel_proj);

        currDetectedOnsetL1  = Tracked_pos_L1.Time(currIdxsL1.start);
        currDetectedOffsetL1 = Tracked_pos_L1.Time(currIdxsL1.end);

        currDetectedOnsetL1  = currDetectedOnsetL1 - config.Speed.timeOffsetForDetectedTemporalWindow;
        if(currDetectedOffsetL1 + config.Speed.timeOffsetForDetectedTemporalWindow > Tracked_pos_L1.Time(end))
            currDetectedOffsetL1 = Tracked_pos_L1.Time(end);
        else
            currDetectedOffsetL1 = currDetectedOffsetL1 + config.Speed.timeOffsetForDetectedTemporalWindow;
        end
        
        % Applying speed cutoff filter 
        Tracked_pos_L1.Filtered_Vel_proj = Tracked_pos_L1.Time >= currDetectedOnsetL1 - 0.0001 & Tracked_pos_L1.Time <= currDetectedOffsetL1 + 0.0001;

        %% Reconstructing l1 length from reconstructed velocity
        currl1real                  = norm(Cone_pos(2,[1 3]) - Cone_pos(1,[1 3]));
        currL1recsmoothed           = trapz(Tracked_pos_L1.Time,Tracked_pos_L1.Smoothed_Vel_proj);
        currL1recsmoothedfiltered   = trapz(Tracked_pos_L1.Time(Tracked_pos_L1.Filtered_Vel_proj),Tracked_pos_L1.Smoothed_Vel_proj(Tracked_pos_L1.Filtered_Vel_proj));
        currL1recUnsmoothed         = trapz(Tracked_pos_L1.Time,Tracked_pos_L1.Vel_proj);
        currL1recUnsmoothedfiltered = trapz(Tracked_pos_L1.Time(Tracked_pos_L1.Filtered_Vel_proj),Tracked_pos_L1.Vel_proj(Tracked_pos_L1.Filtered_Vel_proj));

        if(currl1real < config.Speed.tresholdForBadParticipantL1Recontruction)
            Group.BadPptIdxs = [Group.BadPptIdxs, pId];
        else
            RecontructedQuantities(1,1) = currl1real;
            RecontructedQuantities(1,2) = currL1recsmoothed;
            RecontructedQuantities(1,3) = currL1recsmoothedfiltered;
            RecontructedQuantities(1,4) = currL1recUnsmoothed;
            RecontructedQuantities(1,5) = currL1recUnsmoothedfiltered;
        end

        maxL1 = max(Tracked_pos_L1.Smoothed_Vel_proj);
        idxMaxL1 = find(Tracked_pos_L1.Smoothed_Vel_proj == maxL1);
        Tracked_pos_L1.ShiftedTime = Tracked_pos_L1.Time - Tracked_pos_L1.Time(idxMaxL1);

        %% Calculating a moving average velocity while walking on l2 segment
        l2Size = height(Tracked_pos_L2);
        Tracked_pos_L2.Vel_X = zeros(l2Size,1);
        Tracked_pos_L2.Vel_Z = zeros(l2Size,1);

        for iL = 2:l2Size

            Tracked_pos_L2.Vel_X(iL) = config.Speed.alpha * Tracked_pos_L2.Vel_X(iL-1) +...
                (1 - config.Speed.alpha) * (Tracked_pos_L2.Pos_X(iL) - Tracked_pos_L2.Pos_X(iL-1)) / (Tracked_pos_L2.Time(iL) - Tracked_pos_L2.Time(iL-1));

            Tracked_pos_L2.Vel_Z(iL) = config.Speed.alpha * Tracked_pos_L2.Vel_Z(iL-1) +...
                (1 - config.Speed.alpha) * (Tracked_pos_L2.Pos_Z(iL) - Tracked_pos_L2.Pos_Z(iL-1)) / (Tracked_pos_L2.Time(iL) - Tracked_pos_L2.Time(iL-1));

        end

        % Calculating projection of calculated velocity over ideal walking
        % direction
        Tracked_pos_L2.Vel_proj = dot([Tracked_pos_L2.Vel_X Tracked_pos_L2.Vel_Z],repmat(s_dir_L2,l2Size,1),2);
        % Gaussian filter to smooth data
        Tracked_pos_L2.Smoothed_Vel_proj = smoothdata(Tracked_pos_L2.Vel_proj,'gaussian',config.Speed.smoothWindow);

        % Using a speed threshold cutoff we can isolate from our
        % reconstruction only the data where the participant actually moved 
        filteredVel_proj = Tracked_pos_L2.Smoothed_Vel_proj > config.Speed.velocityCutoff;
        [~,currIdxsL2] = findWalkingSegment(filteredVel_proj);

        currDetectedOnsetL2  = Tracked_pos_L2.Time(currIdxsL2.start);
        currDetectedOffsetL2 = Tracked_pos_L2.Time(currIdxsL2.end);

        currDetectedOnsetL2  = currDetectedOnsetL2 - config.Speed.timeOffsetForDetectedTemporalWindow;
        if(currDetectedOffsetL2 + config.Speed.timeOffsetForDetectedTemporalWindow > Tracked_pos_L2.Time(end))
            currDetectedOffsetL2 = Tracked_pos_L2.Time(end);
        else
            currDetectedOffsetL2 = currDetectedOffsetL2 + config.Speed.timeOffsetForDetectedTemporalWindow;
        end

        % Applying speed cutoff filter 
        Tracked_pos_L2.Filtered_Vel_proj = Tracked_pos_L2.Time >= currDetectedOnsetL2 - 0.0001 & Tracked_pos_L2.Time <= currDetectedOffsetL2 + 0.0001;

        %% Reconstructing l2 length from reconstructed velocity
        currl2real = norm(Cone_pos(3,[1 3]) - Cone_pos(2,[1 3]));
        currL2recsmoothed           = trapz(Tracked_pos_L2.Time,Tracked_pos_L2.Smoothed_Vel_proj);
        currL2recsmoothedfiltered   = trapz(Tracked_pos_L2.Time(Tracked_pos_L2.Filtered_Vel_proj),Tracked_pos_L2.Smoothed_Vel_proj(Tracked_pos_L2.Filtered_Vel_proj));
        currL2recUnsmoothed         = trapz(Tracked_pos_L2.Time,Tracked_pos_L2.Vel_proj);
        currL2recUnsmoothedfiltered = trapz(Tracked_pos_L2.Time(Tracked_pos_L2.Filtered_Vel_proj),Tracked_pos_L2.Vel_proj(Tracked_pos_L2.Filtered_Vel_proj));

        RecontructedQuantities(1,6)  = currl2real;
        RecontructedQuantities(1,7)  = currL2recsmoothed;
        RecontructedQuantities(1,8)  = currL2recsmoothedfiltered;
        RecontructedQuantities(1,9)  = currL2recUnsmoothed;
        RecontructedQuantities(1,10) = currL2recUnsmoothedfiltered;

        maxL2 = max(Tracked_pos_L2.Smoothed_Vel_proj);
        idxMaxL2 = find(Tracked_pos_L2.Smoothed_Vel_proj == maxL2);
        Tracked_pos_L2.ShiftedTime = Tracked_pos_L2.Time - Tracked_pos_L2.Time(idxMaxL2);

        % Calculating times from reconstructed tracking data
        RecontructedQuantities(1,11) = max([currDetectedOffsetL1 - currDetectedOnsetL1,0]); % time spent on l1
        RecontructedQuantities(1,12) = max([currDetectedOnsetL2 - currDetectedOffsetL1, 0]); % time spent standing at cone 2
        RecontructedQuantities(1,13) = max([currDetectedOffsetL2 - currDetectedOnsetL2, 0]); % time spent on l2

        Group.TrackedL1{1,pId}{trialId,1} = Tracked_pos_L1;
        Group.TrackedL2{1,pId}{trialId,1} = Tracked_pos_L2;
        participantRecontructedQuantities = [participantRecontructedQuantities; RecontructedQuantities];

    end

    % Saving output
    Group.Reconstructed{1,pId} = array2table(participantRecontructedQuantities,"VariableNames",...
        {'L1Real' 'L1Smoothed' 'L1SmoothedFiltered' 'L1Unsmoothed' 'L1UnsmoothedFiltered'...
        'L2Real' 'L2Smoothed' 'L2SmoothedFiltered' 'L2Unsmoothed' 'L2UnsmoothedFiltered'...
        'T_L1' 'T_Standing' 'T_L2'});

    Group.Reconstructed{1,pId}.InboundBodyRotation = InboundAngularRotation;
    Group.Reconstructed{1,pId}.SignBodyRotation = SignInboundBodyRotation;
    Group.Reconstructed{1,pId}.InferredReturnAngle = InferredReturnAngle;
    Group.Reconstructed{1,pId}.SignInferredReturnAngle = SignInferredReturnAngle;
    Group.Reconstructed{1,pId}.OobReturnAngle = OobReturnAngle;
    Group.Reconstructed{1,pId}.SignOobReturnAngle = SignOobReturnAngle;
    Group.Reconstructed{1,pId}.RealReturnAngle = RealReturnAngle;
end

Group.overRotations = array2table(Group.overRotations,"VariableNames", {'ParticipantID' 'TrialID'});
Group.BadPptIdxs = unique(Group.BadPptIdxs);

outGroup = Group;

    function [y, idxs] = findWalkingSegment(x)
    % Find the biggest train of 1s within a logical array   
        
        demA = 0;
        demB = 0;
        a = zeros(length(x),1);
        b = zeros(length(x),1);
        y = x;
        for i=length(x):-1:1
            if x(i) == 1
                demA = demA + 1;
            else
                demA = 0;
            end
    
            if(x(length(x) - i + 1) == 1)
                demB = demB + 1;
            else
                demB = 0;
            end
    
            a(i) = demA;
            b(length(x) - i + 1) = demB;
        end
    
        xStart = find(a == max(a),1);
        xEnd = find(b == max(b),1);
    
        %Cutting before start and after end
        y(1:xStart-1) = 0;
        y(xEnd+1:length(x)) = 0;
    
        idxs.start = xStart;
        idxs.end = xEnd;
    
    end

end

