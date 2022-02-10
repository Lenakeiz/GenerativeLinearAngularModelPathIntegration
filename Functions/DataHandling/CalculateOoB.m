function [ Errors,  OoBInfo] = CalculateOoB(FlagPos, TrigPos, OutOfBoundPos, CondTable)
        
    debug_plotting = 0;
    
    %% Getting cone positions
    flagPosMat = cell2mat(FlagPos);
    firstConePositions = flagPosMat(1 : 3 : end, [1,3]);
    secondConePositions = flagPosMat(2 : 3 : end, [1,3]);
    thirdConePositions = flagPosMat(3 : 3 : end, [1,3]);

    %% Getting triggered position excluding y values
    triggerPositions = cell2mat(TrigPos);
    triggerPositions = triggerPositions(:,[1,3]);

    %% Getting triggered position excluding y values
    outofboundpos = cell2mat(OutOfBoundPos);
    outofboundpos = outofboundpos(:,[1,3]);

    %% Calculating absolute error distance as sqrt((x2-x1)^2 + (z2-z1)^2)
    AbsErrorDistance = sqrt(sum((triggerPositions - firstConePositions).^2,2));
    AbsWalkedSegment = sqrt(sum((triggerPositions - thirdConePositions).^2,2));
    OoBLength = zeros(size(triggerPositions,1),1);

    OoBInfo.ReconstructedOoB = OutOfBoundPos;
    OoBInfo.IsRecontructedFromOverallMean = zeros(size(triggerPositions,1),1);

    meanReturnSegment = nanmean(AbsWalkedSegment(CondTable.OutOfBound == 0));
    meanReturnSegmentNoChange = nanmean(AbsWalkedSegment(CondTable.Condition == 1 & CondTable.OutOfBound == 0));
    meanReturnSegmentNoDistalCues = nanmean(AbsWalkedSegment(CondTable.Condition == 2 & CondTable.OutOfBound == 0));
    meanReturnSegmentNoOpticFlow = nanmean(AbsWalkedSegment(CondTable.Condition == 3 & CondTable.OutOfBound == 0));

    idx = find(CondTable.OutOfBound == 1);
    


    for i=1:size(idx,1)

        disp("Plotting trial " + idx(i));
        currOoB = outofboundpos(idx(i),:);
        currOoBLength = sqrt(sum((currOoB - thirdConePositions(idx(i),:)).^2,2));
        OoBLength(idx(i)) = currOoBLength;
        directionOoB = currOoB - thirdConePositions(idx(i),:);
        directionOoB = directionOoB./norm(directionOoB);

        inferredOoBLength = nan;

        switch CondTable.Condition(idx(i))
            case 1
                inferredOoBLength = meanReturnSegmentNoChange;
            case 2
                inferredOoBLength = meanReturnSegmentNoDistalCues;
            case 3
                inferredOoBLength = meanReturnSegmentNoOpticFlow;
        end

        if(isnan(inferredOoBLength))
            %Flagging the fact we are using the overall return mean for
            %this calculation (it means the participant did not complete
            %any successfull trial for that condition
            inferredOoBLength = meanReturnSegment;
            OoBInfo.IsRecontructedFromOverallMean(idx(i),1) = 1;
            disp("Reconstruction with overall mean " + idx(i));
        end 

        recontructedVert = thirdConePositions(idx(i),:) + directionOoB*inferredOoBLength;
 
%         if(currOoBLength < inferredOoBLength)
%             OoBInfo.ReconstructedOoB{idx(i),:} = [recontructedVert(1) 0 recontructedVert(2)];
%         else
%             disp("Discarding current inferred location");
%             OoBInfo.ReconstructedOoB{idx(i),:} = [currOoB(1) 0 currOoB(2)];
%         end
        
        %for all the OoB trials, replace the leg length with the infered
        %one
        OoBInfo.ReconstructedOoB{idx(i),:} = [recontructedVert(1) 0 recontructedVert(2)];

        % Debug
        if(debug_plotting == 1)
            a = [firstConePositions(idx(i),:);secondConePositions(idx(i),:);thirdConePositions(idx(i),:);currOoB;OoBInfo.ReconstructedOoB{idx(i)}(:,[1 3])]; scatter(a(:,1), a(:,2)); hold on; plot(a(:,1), a(:,2));
            ts = ["X_0";"X_1";"X_2";"OoB";"OoB_i"]; for j = 1:5; text(a(j,1),a(j,2),ts(j,1)); end; hold off; xlim([-1 7]); ylim([-4 4]); axis square;
            titletype = ["No Change"; "No Distal Cues"; "No Optic Flow"]; title(titletype(CondTable.Condition(idx(i)))); drawnow; waitforbuttonpress;            
        end
    end
    
    Errors.AbsErrorDistance      = AbsErrorDistance;
    Errors.AbsWalkedSegment      = AbsWalkedSegment;
    Errors.OoBLength             = OoBLength;
    Errors.MeanLastWalkedSegment = array2table([meanReturnSegment meanReturnSegmentNoChange meanReturnSegmentNoDistalCues meanReturnSegmentNoOpticFlow],'VariableNames',{'Total', 'NoChange', 'NoDistalCue', 'NoOpticFlow'});
end