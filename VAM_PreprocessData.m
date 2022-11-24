% %% Preprocessing Young Data
% disp("%%%%%%%%%%%%%%% Preprocessing young controls %%%%%%%%%%%%%%%");
% config.Speed.tresholdForBadParticipantL1Recontruction   = 1.55;    % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
% YoungControls = TransformPaths(YoungControls);%transform data
% YoungControls   = CalculateTrackingPath(YoungControls, config);
% ManuallyScoringYoung;
%% Preprocessing Older Control Data
disp("%%%%%%%%%%%%%%% Preprocessing healthy controls %%%%%%%%%%%%%%%");
config.Speed.tresholdForBadParticipantL1Recontruction   = 2.0;    % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
HealthyControls = TransformPaths(HealthyControls);%transform data
HealthyControls   = CalculateTrackingPath(HealthyControls, config);
ManuallyScoringHealthyOld;
% %% Preprocessing MCIPos Control Data
% disp("%%%%%%%%%%%%%%% Starting fit for mci positive %%%%%%%%%%%%%%%");
% config.Speed.tresholdForBadParticipantL1Recontruction   = 0.0;    % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
% MCIPos = TransformPaths(MCIPos);%transform data
% MCIPos   = CalculateTrackingPath(MCIPos, config);
% ManuallyScoringMCIPos;
% %% Preprocessing MCI Neg
% disp("%%%%%%%%%%%%%%% Preprocessing mci negative %%%%%%%%%%%%%%%");
% config.Speed.tresholdForBadParticipantL1Recontruction   = 0.0; 
% MCINeg = TransformPaths(MCINeg);%transform data
% MCINeg   = CalculateTrackingPath(MCINeg, config);
% ManuallyScoringMCINeg;
% %%
% disp("%%%%%%%%%%%%%%% Preprocessing mci unknown %%%%%%%%%%%%%%%");
% config.Speed.tresholdForBadParticipantL1Recontruction   = 0.0; 
% MCIUnk = TransformPaths(Unknown);%transform data
% MCIUnk   = CalculateTrackingPath(MCIUnk, config);
% ManuallyScoringMCIUnk;
% clear Unknown % removing duplicate
% 
% disp("%%%%%%%%%%%%%%% Preprocessing complete %%%%%%%%%%%%%%%");
