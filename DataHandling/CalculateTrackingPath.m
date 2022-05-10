function [outGroup] = CalculateTrackingPath(Group, config)
    %% Plotting some example data
    
    pSize = size(Group.Path,2);
    
    Group.TrackedL1 = {};
    Group.TrackedL2 = {};
    Group.Reconstructed = {};
    Group.BadPptIdxs = [];

    % For angular calculation. Calculating the signed angle between two
    % vectors (please refer to
    % https://wumbo.net/formula/angle-between-two-vectors-2d/)
    anglebetween = @(va,vb) atan2d(va(:,1).*vb(:,2) - va(:,2).*vb(:,1), va(:,1).*vb(:,1) + va(:,2).*vb(:,2));

    % Delta t multiplier is used to integrate every dtMultiplier delta
    % times instead of 100 ms (in this second every 500 ms)
    dtMultiplier = config.TrackedInboundAngularDeltaT;

    for pId = 1:pSize
        
        trialSize = size(Group.Path{1,pId},1);

        Group.TrackedL1{1,pId} = {};
        Group.TrackedL2{1,pId} = {};
   
        InboundAngularRotation = [];

        participantRecontructedQuantities = [];

        % Calculate the mean time between making cone 1 disappear and spawning cone
        % two (same between 2 and 3)
        % Spoiler alert, is always 2 seconds
        DiffTrigSpawn12 = cell2mat(Group.FlagSpawnTime{1,pId}(:,2)) - cell2mat(Group.FlagTrigTimes{1,pId}(:,1));
        DiffTrigSpawn23 = cell2mat(Group.FlagSpawnTime{1,pId}(:,3)) - cell2mat(Group.FlagTrigTimes{1,pId}(:,2));
    
        for trialId = 1:trialSize
            
            Cone_pos    = Group.FlagPos{1,pId}{trialId,1};
            Tracked_pos = Group.Path{1,pId}{trialId,1};
            Tracked_pos = array2table(Tracked_pos,"VariableNames",{'Time' 'Pos_X' 'Pos_Y' 'Pos_Z' 'Forward_X' 'Forward_Y' 'Forward_Z'});
            
            %% Calculating angular rotations from tracking data
            % For angular calculations getting only the tracked position after reaching cone 3
            Angular_pos = Tracked_pos(Tracked_pos.Time > Group.FlagTrigTimes{1,pId}{trialId,3},:);
            % Getting only the portion of angular data we are interested in
            Angular_pos = removevars(Tracked_pos,{'Pos_X' 'Pos_Y' 'Pos_Z' 'Forward_Y'});
            % Getting a reading only after dtMultiplier * dt 
            Angular_pos = Tracked_pos(1:dtMultiplier:end,:);
            % Calculating the turning angle for this trial by summing up
            % all of the angular dts
            tracking_size = height(Angular_pos);
            ang_rotation = 0;
            for trackId = 2:tracking_size
                ang_rotation = ang_rotation + anglebetween([Angular_pos.Forward_X(trackId) Angular_pos.Forward_Z(trackId)],[Angular_pos.Forward_X(trackId-1) Angular_pos.Forward_Z(trackId-1)]);
            end
            InboundAngularRotation = [InboundAngularRotation;ang_rotation];

            %% Calculating speed and tracked distances
            % VARIABLE NAMES --> l1real l1realsmoothed
            % l1realsmoothedfilterd l1realunsmoothed
            % l1realunsmoothedfiltered l2real l2realsmoothed
            % l2realsmoothedfilterd l2realunsmoothed
            % l2realunsmoothedfiltered t1(between l1 onset - offset)
            % t2(between l1 offset - l2 onset) t3(between l2 onset -
            % offset)
            RecontructedQuantities = nan(1,13);

            % Filtering positions
            % Excluding tracking before reaching cone 1
            Tracked_pos    = Tracked_pos(Tracked_pos.Time > Group.FlagTrigTimes{1,pId}{trialId,1},:);
            
            % Extracting L1 walking path (between flag spawn time 2 and reaching cone 2)

            Tracked_pos_L1 = Tracked_pos(Tracked_pos.Time > Group.FlagSpawnTime{1,pId}{trialId,2} &...
                                        Tracked_pos.Time < Group.FlagTrigTimes{1,pId}{trialId,2} + mean(DiffTrigSpawn12) + config.Speed.timeOffsetAfterFlagReach, :);
%             Tracked_pos_L1 = Tracked_pos(Tracked_pos.Time > Group.FlagTrigTimes{1,pId}{trialId,1} &...
%                             Tracked_pos.Time < Group.FlagTrigTimes{1,pId}{trialId,2} + config.Speed.TOffsetAfterFlagReach, :);
            
            % Extracting L1 walking path (between flag spawn time 3 and reaching cone 3)
            Tracked_pos_L2 = Tracked_pos(Tracked_pos.Time > Group.FlagSpawnTime{1,pId}{trialId,3} &...
                                         Tracked_pos.Time < Group.FlagTrigTimes{1,pId}{trialId,3} + mean(DiffTrigSpawn23) + config.Speed.timeOffsetAfterFlagReach, :);

%             Tracked_pos_L2 = Tracked_pos(Tracked_pos.Time > Group.FlagTrigTimes{1,pId}{trialId,2} &...
%                                          Tracked_pos.Time < Group.FlagTrigTimes{1,pId}{trialId,3} + config.Speed.TOffsetAfterFlagReach, :);

            % Calculating ideal direction of walking
            s_dir_L1 = (Cone_pos(2,[1 3]) - Cone_pos(1,[1 3])) / norm(Cone_pos(2,[1 3]) - Cone_pos(1,[1 3])); 
            s_dir_L2 = (Cone_pos(3,[1 3]) - Cone_pos(2,[1 3])) / norm(Cone_pos(3,[1 3]) - Cone_pos(2,[1 3])); 

            l1Size = height(Tracked_pos_L1);
            Tracked_pos_L1.Vel_X = zeros(l1Size,1);
            Tracked_pos_L1.Vel_Z = zeros(l1Size,1);

            for iL = 2:l1Size
            
                Tracked_pos_L1.Vel_X(iL) = config.Speed.alpha * Tracked_pos_L1.Vel_X(iL-1) +...
                                            (1 - config.Speed.alpha) * (Tracked_pos_L1.Pos_X(iL) - Tracked_pos_L1.Pos_X(iL-1)) / (Tracked_pos_L1.Time(iL) - Tracked_pos_L1.Time(iL-1)); 
            
                Tracked_pos_L1.Vel_Z(iL) = config.Speed.alpha * Tracked_pos_L1.Vel_Z(iL-1) +...
                                            (1 - config.Speed.alpha) * (Tracked_pos_L1.Pos_Z(iL) - Tracked_pos_L1.Pos_Z(iL-1)) / (Tracked_pos_L1.Time(iL) - Tracked_pos_L1.Time(iL-1));  
            
            end

            % Calculating projecting over the cone2-1 direction
            Tracked_pos_L1.Vel_proj = dot([Tracked_pos_L1.Vel_X Tracked_pos_L1.Vel_Z],repmat(s_dir_L1,l1Size,1),2);
            % Using a gaussian filter to smooth even more the data
            Tracked_pos_L1.Smoothed_Vel_proj = smoothdata(Tracked_pos_L1.Vel_proj,'gaussian',config.Speed.smoothWindow);

            % Calculating velocity cutoff
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

            % The additional 0.0001 is to make sure that the resolution
            % gets the initial detected value (comparing to eps is still
            % not enough).
            Tracked_pos_L1.Filtered_Vel_proj = Tracked_pos_L1.Time >= currDetectedOnsetL1 - 0.0001 & Tracked_pos_L1.Time <= currDetectedOffsetL1 + 0.0001;

            % Calculating the reconstructed quantities for L1
            currl1real                  = norm(Cone_pos(2,[1 3]) - Cone_pos(1,[1 3]));
            currL1recsmoothed           = trapz(Tracked_pos_L1.Time,Tracked_pos_L1.Smoothed_Vel_proj);
            currL1recsmoothedfiltered   = trapz(Tracked_pos_L1.Time(Tracked_pos_L1.Filtered_Vel_proj),Tracked_pos_L1.Smoothed_Vel_proj(Tracked_pos_L1.Filtered_Vel_proj));
            currL1recUnsmoothed         = trapz(Tracked_pos_L1.Time,Tracked_pos_L1.Vel_proj);
            currL1recUnsmoothedfiltered = trapz(Tracked_pos_L1.Time(Tracked_pos_L1.Filtered_Vel_proj),Tracked_pos_L1.Vel_proj(Tracked_pos_L1.Filtered_Vel_proj));

            % This happens only for L1 at the moment
            if(currl1real < config.Speed.tresholdForBadParticipantL1Recontruction)
                Group.BadPptIdxs = [Group.BadPptIdxs, pId];
            else
                RecontructedQuantities(1,1) = currl1real;
                RecontructedQuantities(1,2) = currL1recsmoothed;
                RecontructedQuantities(1,3) = currL1recsmoothedfiltered;
                RecontructedQuantities(1,4) = currL1recUnsmoothed;
                RecontructedQuantities(1,5) = currL1recUnsmoothedfiltered;
            end

            %Calculating the shifted position for the smoothed velocity
            %projection
            maxL1 = max(Tracked_pos_L1.Smoothed_Vel_proj);
            idxMaxL1 = find(Tracked_pos_L1.Smoothed_Vel_proj == maxL1);
            Tracked_pos_L1.ShiftedTime = Tracked_pos_L1.Time - Tracked_pos_L1.Time(idxMaxL1);

            l2Size = height(Tracked_pos_L2);
            Tracked_pos_L2.Vel_X = zeros(l2Size,1);
            Tracked_pos_L2.Vel_Z = zeros(l2Size,1);

            for iL = 2:l2Size
            
                Tracked_pos_L2.Vel_X(iL) = config.Speed.alpha * Tracked_pos_L2.Vel_X(iL-1) +...
                                    (1 - config.Speed.alpha) * (Tracked_pos_L2.Pos_X(iL) - Tracked_pos_L2.Pos_X(iL-1)) / (Tracked_pos_L2.Time(iL) - Tracked_pos_L2.Time(iL-1)); 
                
                Tracked_pos_L2.Vel_Z(iL) = config.Speed.alpha * Tracked_pos_L2.Vel_Z(iL-1) +...
                                (1 - config.Speed.alpha) * (Tracked_pos_L2.Pos_Z(iL) - Tracked_pos_L2.Pos_Z(iL-1)) / (Tracked_pos_L2.Time(iL) - Tracked_pos_L2.Time(iL-1));  
            
            end

            % Calculating projecting over the cone3-1 direction
            Tracked_pos_L2.Vel_proj = dot([Tracked_pos_L2.Vel_X Tracked_pos_L2.Vel_Z],repmat(s_dir_L2,l2Size,1),2);
            Tracked_pos_L2.Smoothed_Vel_proj = smoothdata(Tracked_pos_L2.Vel_proj,'gaussian',config.Speed.smoothWindow);

            % Calculating velocity cutoff
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

            Tracked_pos_L2.Filtered_Vel_proj = Tracked_pos_L2.Time >= currDetectedOnsetL2 - 0.0001 & Tracked_pos_L2.Time <= currDetectedOffsetL2 + 0.0001;

            % Calculating the reconstructed quantities for L2
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

            % Calculating times
            RecontructedQuantities(1,11) = max([currDetectedOffsetL1 - currDetectedOnsetL1,0]);
            RecontructedQuantities(1,12) = max([currDetectedOnsetL2 - currDetectedOffsetL1, 0]);
            RecontructedQuantities(1,13) = max([currDetectedOffsetL2 - currDetectedOnsetL2, 0]);

            %Calculating the shifted position for the smoothed velocity
            %projection
            maxL2 = max(Tracked_pos_L2.Smoothed_Vel_proj);
            idxMaxL2 = find(Tracked_pos_L2.Smoothed_Vel_proj == maxL2);
            Tracked_pos_L2.ShiftedTime = Tracked_pos_L2.Time - Tracked_pos_L2.Time(idxMaxL2);

            Group.TrackedL1{1,pId}{trialId,1} = Tracked_pos_L1;
            Group.TrackedL2{1,pId}{trialId,1} = Tracked_pos_L2;
            participantRecontructedQuantities = [participantRecontructedQuantities; RecontructedQuantities];

        end

        % VARIABLE NAMES --> l1real l1realsmoothed
        % l1realsmoothedfilterd l1realunsmoothed
        % l1realunsmoothedfiltered l2real l2realsmoothed
        % l2realsmoothedfilterd l2realunsmoothed
        % l2realunsmoothedfiltered t1(between l1 onset - offset)
        % t2(between l1 offset - l2 onset) t3(between l2 onset -
        % offset)
        Group.Reconstructed{1,pId} = array2table(participantRecontructedQuantities,"VariableNames",...
            {'L1Real' 'L1Smoothed' 'L1SmoothedFiltered' 'L1Unsmoothed' 'L1UnsmoothedFiltered'...
             'L2Real' 'L2Smoothed' 'L2SmoothedFiltered' 'L2Unsmoothed' 'L2UnsmoothedFiltered'...
             'T_L1' 'T_Standing' 'T_L2'});
        % Adding angular calculation
        Group.Reconstructed{1,pId}.InboundRotation = InboundAngularRotation;
    
    end

    Group.BadPptIdxs = unique(Group.BadPptIdxs);

    outGroup = Group;

    function [y, idxs] = findWalkingSegment(x)
        
        % This function returns the finds the index for the maximum maximum i
        demA = 0;
        demB = 0;
        a = zeros(length(x),1);
        b = zeros(length(x),1);
        y = x;
            for i=length(x):-1:1
                %Traversing from the back
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

