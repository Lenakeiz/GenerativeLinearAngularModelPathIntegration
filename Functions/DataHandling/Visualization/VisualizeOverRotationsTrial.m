function VisualizeOverRotationsTrial(Group, reproductionSpeed, groupName, videoplot)
% Andrea Castegnaro, UCL, 2022 uceeaca@ucl.ac.uk
% Helper function to visualize in realtime all of the path that have been
% indicated to show an over rotation of the participant, based on the
% tracked body position. Over rotations occurred when at cone 3 the
% participant turned the body on the long side of the rotation rather then
% taking the shorter angle. Additional information can be found in the
% CalculateTrackingPath function.
% ===================================================================================

if(~isfield(Group,'overRotations'))
    disp('Please run CalculateTrckingPath first');
    return;
end

for idx = 1:height(Group.overRotations)

    VisualizeRealtimeTrackingData(Group,Group.overRotations.ParticipantID(idx), Group.overRotations.TrialID(idx), reproductionSpeed, 'cutconethree', true, 'groupname', groupName, "videoplot", videoplot);

end

end