function [ Errors,  OoBInfo] = CalculateErrors(FlagPos, TrigPos, OutOfBoundPos, CondTable)
%% CalculateErrors 
% Andrea Castegnaro, UCL, 2022 uceeaca@ucl.ac.uk
% Calculate absolute distance error for each set of trials divided by
% environmental condition. If out of bound (OoB) trial calculates the
% return distance set as the mean value of all good trials for that
% environmental condition
% ===================================================================================

% local variable for quick visualization. 
debug_plotting = 0;

% Getting cone positions
flagPosMat = cell2mat(FlagPos);
firstConePositions = flagPosMat(1 : 3 : end, [1,3]);
secondConePositions = flagPosMat(2 : 3 : end, [1,3]);
thirdConePositions = flagPosMat(3 : 3 : end, [1,3]);

% Getting triggered position excluding y values
triggerPositions = cell2mat(TrigPos);
triggerPositions = triggerPositions(:,[1,3]);

% Getting OoB position excluding y values
outofboundpos = cell2mat(OutOfBoundPos);
outofboundpos = outofboundpos(:,[1,3]);

% Calculating absolute error distance
AbsErrorDistance = sqrt(sum((triggerPositions - firstConePositions).^2,2));
AbsWalkedSegment = sqrt(sum((triggerPositions - thirdConePositions).^2,2));
OoBLength = zeros(size(triggerPositions,1),1);

% Setting output
OoBInfo.ReconstructedOoB = OutOfBoundPos;
OoBInfo.IsRecontructedFromOverallMean = zeros(size(triggerPositions,1),1);

meanReturnSegment = mean(AbsWalkedSegment(CondTable.OutOfBound == 0),"omitnan");
meanReturnSegmentNoChange = mean(AbsWalkedSegment(CondTable.Condition == 1 & CondTable.OutOfBound == 0),"omitnan");
meanReturnSegmentNoDistalCues = mean(AbsWalkedSegment(CondTable.Condition == 2 & CondTable.OutOfBound == 0),"omitnan");
meanReturnSegmentNoOpticFlow = mean(AbsWalkedSegment(CondTable.Condition == 3 & CondTable.OutOfBound == 0),"omitnan");

idx = find(CondTable.OutOfBound == 1);

for i=1:size(idx,1)

    currOoB = outofboundpos(idx(i),:);
    currOoBLength = sqrt(sum((currOoB - thirdConePositions(idx(i),:)).^2,2));
    OoBLength(idx(i)) = currOoBLength;
    directionOoB = currOoB - thirdConePositions(idx(i),:);
    directionOoB = directionOoB./norm(directionOoB);

    % For OoB trials we set the distance error as the mean distance for all
    % of the trials for that conditions that are valid
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
        % If the participant has no valid trials for that condition we are
        % going to the use the mean absolute distance error on the return
        % path from all of the trials. We flag this condition for later.
        inferredOoBLength = meanReturnSegment;
        OoBInfo.IsRecontructedFromOverallMean(idx(i),1) = 1;
        disp("Reconstruction with overall mean " + idx(i));
    end

    recontructedVert = thirdConePositions(idx(i),:) + directionOoB*inferredOoBLength;

    % For all the OoB trials, replace the leg length with the inferred
    % One
    OoBInfo.ReconstructedOoB{idx(i),:} = [recontructedVert(1) 0 recontructedVert(2)];

    % Quick debug plot. Will stop execution code. Variable is set at the
    % beginning of the script.
    if(debug_plotting == 1)
        disp("Plotting trial " + idx(i));
        a = [firstConePositions(idx(i),:);secondConePositions(idx(i),:);thirdConePositions(idx(i),:);currOoB;OoBInfo.ReconstructedOoB{idx(i)}(:,[1 3])]; scatter(a(:,1), a(:,2)); hold on; plot(a(:,1), a(:,2));
        ts = ["X_0";"X_1";"X_2";"OoB";"OoB_i"]; for j = 1:5; text(a(j,1),a(j,2),ts(j,1)); end; hold off; xlim([-1 7]); ylim([-4 4]); axis square;
        titletype = ["No Change"; "No Distal Cues"; "No Optic Flow"]; title(titletype(CondTable.Condition(idx(i)))); drawnow; waitforbuttonpress;
    end
end

% Setting output
Errors.AbsErrorDistance      = AbsErrorDistance;
Errors.AbsWalkedSegment      = AbsWalkedSegment;
Errors.OoBLength             = OoBLength;
Errors.MeanLastWalkedSegment = array2table([meanReturnSegment meanReturnSegmentNoChange meanReturnSegmentNoDistalCues meanReturnSegmentNoOpticFlow],'VariableNames',{'Total', 'NoChange', 'NoDistalCue', 'NoOpticFlow'});
end