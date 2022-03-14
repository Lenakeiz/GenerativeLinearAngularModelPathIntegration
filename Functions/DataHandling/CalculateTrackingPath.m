function [outGroup] = CalculateTrackingPath(Group, config)
    %% Plotting some example data
    
    pSize = size(Group.Path,2);
    
    Group.TrackedL1 = {};
    Group.TrackedL2 = {};
    
    for pId = 1:pSize
        
        Group.TrackedL1{1,pId} = {};
        Group.TrackedL2{1,pId} = {};

        trialSize = size(Group.Path{1,pId},1);
    
        % Calculate the mean time between making cone 1 disappear and spawning cone
        % two (same between 2 and 3)
        % Spoiler alert, is always 2 seconds
        DiffTrigSpawn12 = cell2mat(Group.FlagSpawnTime{1,pId}(:,2)) - cell2mat(Group.FlagTrigTimes{1,pId}(:,1));
        DiffTrigSpawn23 = cell2mat(Group.FlagSpawnTime{1,pId}(:,3)) - cell2mat(Group.FlagTrigTimes{1,pId}(:,2));
    
        for trialId = 1:trialSize
            
            Cone_pos    = Group.FlagPos{1,pId}{trialId,1};
            Tracked_pos = Group.Path{1,pId}{trialId,1};
            Tracked_pos = array2table(Tracked_pos,"VariableNames",{'Time' 'Pos_X' 'Pos_Y' 'Pos_Z' 'Euler_X' 'Euler_Y' 'Euler_Z'});
            
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
            Tracked_pos_L1.Smoothed_Vel_proj = smoothdata(Tracked_pos_L1.Vel_proj,'gaussian',config.Speed.smoothWindow);

            % Calculating velocity cutoff
            filteredVel_proj = Tracked_pos_L1.Smoothed_Vel_proj > config.Speed.velocityCutoff;
            [~,currIdxs] = findWalkingSegment(filteredVel_proj);

            currDetectedOnset  = Tracked_pos_L1.Time(currIdxs.start);
            currDetectedOffset = Tracked_pos_L1.Time(currIdxs.end);

            currDetectedOnset  = currDetectedOnset - config.Speed.timeOffsetForDetectedTemporalWindow;
            if(currDetectedOffset + config.Speed.timeOffsetForDetectedTemporalWindow > Tracked_pos_L1.Time(end))
                currDetectedOffset = Tracked_pos_L1.Time(end);
            else
                currDetectedOffset = currDetectedOffset + config.Speed.timeOffsetForDetectedTemporalWindow;
            end

            % The additional 0.0001 is to make sure that the resolution
            % gets the initial detected value (comparing to eps is still
            % not enough).
            Tracked_pos_L1.Filthered_Vel_proj = Tracked_pos_L1.Time >= currDetectedOnset - 0.0001 & Tracked_pos_L1.Time <= currDetectedOffset + 0.0001;

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
            [~,currIdxs] = findWalkingSegment(filteredVel_proj);

            currDetectedOnset  = Tracked_pos_L2.Time(currIdxs.start);
            currDetectedOffset = Tracked_pos_L2.Time(currIdxs.end);

            currDetectedOnset  = currDetectedOnset - config.Speed.timeOffsetForDetectedTemporalWindow;
            if(currDetectedOffset + config.Speed.timeOffsetForDetectedTemporalWindow > Tracked_pos_L2.Time(end))
                currDetectedOffset = Tracked_pos_L2.Time(end);
            else
                currDetectedOffset = currDetectedOffset + config.Speed.timeOffsetForDetectedTemporalWindow;
            end

            Tracked_pos_L2.Filthered_Vel_proj = Tracked_pos_L2.Time >= currDetectedOnset - 0.0001 & Tracked_pos_L2.Time <= currDetectedOffset + 0.0001;

            %Calculating the shifted position for the smoothed velocity
            %projection
            maxL2 = max(Tracked_pos_L2.Smoothed_Vel_proj);
            idxMaxL2 = find(Tracked_pos_L2.Smoothed_Vel_proj == maxL2);
            Tracked_pos_L2.ShiftedTime = Tracked_pos_L2.Time - Tracked_pos_L2.Time(idxMaxL2);

            Group.TrackedL1{1,pId}{trialId,1} = Tracked_pos_L1;
            Group.TrackedL2{1,pId}{trialId,1} = Tracked_pos_L2;

        end
    
    end

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
    
        xStart = find(a == max(a));
        xEnd = find(b == max(b));
    
        %Cutting before start and after end
        y(1:xStart-1) = 0;
        y(xEnd+1:length(x)) = 0;

        idxs.start = xStart;
        idxs.end = xEnd;
    
    end

end

