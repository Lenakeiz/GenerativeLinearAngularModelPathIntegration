function [outGroup] = CalculateTrackingPath(Group)
    %% Plotting some example data
    
    pSize = size(Group.Path,2);
    
    Group.L1_Paths = {};
    Group.L2_Paths = {};
    
    for pId = 1:pSize
        
        Group.L1_Paths{1,pId} = {};
        Group.L2_Paths{1,pId} = {};

        trialSize = size(Group.Path{1,pId},1);
    
        % Calculate the mean time between making cone 1 disappear and spawning cone
        % two (same between 2 and 3)
        DiffTrigSpawn12 = cell2mat(Group.FlagSpawnTime{1,pId}(:,2)) - cell2mat(Group.FlagTrigTimes{1,pId}(:,1));
        DiffTrigSpawn23 = cell2mat(Group.FlagSpawnTime{1,pId}(:,3)) - cell2mat(Group.FlagTrigTimes{1,pId}(:,2));
    
        for trialId = 1:trialSize
            
            Tracked_pos = Group.Path{1,pId}{trialId,1};
            Tracked_pos = array2table(Tracked_pos,"VariableNames",{'Time' 'Pos_X' 'Pos_Y' 'Pos_Z' 'Euler_X' 'Euler_Y' 'Euler_Z'});
            
            % Filtering positions
            % Excluding tracking before reaching cone 1
            Tracked_pos    = Tracked_pos(Tracked_pos.Time > Group.FlagTrigTimes{1,pId}{trialId,1},:);
            
            % Extracting L1 walking path (between flag spawn time 2 and reaching cone 2)
            Tracked_pos_L1 = Tracked_pos(Tracked_pos.Time > Group.FlagSpawnTime{1,pId}{trialId,2} &...
                                         Tracked_pos.Time < Group.FlagTrigTimes{1,pId}{trialId,2} + mean(DiffTrigSpawn12)/2 , :);
            
            % Extracting L1 walking path (between flag spawn time 2 and reaching cone 2)
            Tracked_pos_L2 = Tracked_pos(Tracked_pos.Time > Group.FlagSpawnTime{1,pId}{trialId,3} &...
                                         Tracked_pos.Time < Group.FlagTrigTimes{1,pId}{trialId,3} + mean(DiffTrigSpawn23)/2, :);
        
            Group.L1_Paths{1,pId}{trialId,1} = Tracked_pos_L1;
            Group.L2_Paths{1,pId}{trialId,1} = Tracked_pos_L2;

            % Calculating velocity for each segment
            

        end
    
    end

    outGroup = Group;

end

