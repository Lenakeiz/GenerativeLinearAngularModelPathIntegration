function VisualizeOverRotationsTrial(Group, reproductionSpeed)

if(~isfield(Group,'overRotations'))
    disp('Please run CalculateTrckingPath first');
    return;
end

for idx = 1:height(Group.overRotations)

    VisualizeRealtimeTrackingData(Group,Group.overRotations.ParticipantID(idx), Group.overRotations.TrialID(idx),reproductionSpeed,'cutconethree',true);

end

end