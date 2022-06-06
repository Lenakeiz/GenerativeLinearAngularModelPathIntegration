%Run Calculate TrackingPath first

%% Participant 1
trialsize = height(MCIPos.Reconstructed{1,1});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIPos,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIPos.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
BadExecution(8,1) = 1;
BadExecution(20,1) = 1;
%Save the information in the Reconstructed Structure
MCIPos.Reconstructed{1,1}.BadExecution = BadExecution;

%% Participant 2
trialsize = height(MCIPos.Reconstructed{1,2});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIPos,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIPos.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
BadExecution(3,1) = 1;
BadExecution(6,1) = 1;
BadExecution(10,1) = 1;
BadExecution(12,1) = 1;
BadExecution(14,1) = 1;
BadExecution(25,1) = 1;
BadExecution(27,1) = 1;

MCIPos.Reconstructed{1,2}.RealReturnAngle(5) = MCIPos.Reconstructed{1,2}.RealReturnAngle(5) + 360;
MCIPos.Reconstructed{1,2}.RealReturnAngle(16) = MCIPos.Reconstructed{1,2}.RealReturnAngle(16) + 360;
%Save the information in the Reconstructed Structure
MCIPos.Reconstructed{1,2}.BadExecution = BadExecution;

%% Participant 3
trialsize = height(MCIPos.Reconstructed{1,3});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIPos,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIPos.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
BadExecution(22,1) = 1;
BadExecution(25,1) = 1;
%Save the information in the Reconstructed Structure
MCIPos.Reconstructed{1,3}.BadExecution = BadExecution;

%% Participant 4
trialsize = height(MCIPos.Reconstructed{1,4});
% Prepare structure
BadExecution = zeros(trialsize,1);

MCIPos.Reconstructed{1,4}.RealReturnAngle(8) = MCIPos.Reconstructed{1,4}.RealReturnAngle(8) + 360;

%Save the information in the Reconstructed Structure
MCIPos.Reconstructed{1,4}.BadExecution = BadExecution;

%% Participant 5
trialsize = height(MCIPos.Reconstructed{1,5});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIPos,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIPos.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
BadExecution(3,1) = 1;

MCIPos.Reconstructed{1,5}.RealReturnAngle(2) = MCIPos.Reconstructed{1,5}.RealReturnAngle(2) + 360;
MCIPos.Reconstructed{1,5}.RealReturnAngle(9) = MCIPos.Reconstructed{1,5}.RealReturnAngle(9) - 360;

MCIPos.Reconstructed{1,5}.BadExecution = BadExecution;

%% Participant 6
trialsize = height(MCIPos.Reconstructed{1,6});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIPos,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIPos.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
MCIPos.Reconstructed{1,6}.RealReturnAngle(8) = MCIPos.Reconstructed{1,6}.RealReturnAngle(8) + 360;

%Save the information in the Reconstructed Structure
MCIPos.Reconstructed{1,6}.BadExecution = BadExecution;

%% Participant 7
trialsize = height(MCIPos.Reconstructed{1,7});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(10,1) = 1;

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIPos,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIPos.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
MCIPos.Reconstructed{1,7}.RealReturnAngle(5) = MCIPos.Reconstructed{1,7}.RealReturnAngle(5) + 360;
MCIPos.Reconstructed{1,7}.RealReturnAngle(13) = MCIPos.Reconstructed{1,7}.RealReturnAngle(13) + 360;

%Save the information in the Reconstructed Structure
MCIPos.Reconstructed{1,7}.BadExecution = BadExecution;

%% Participant 8
trialsize = height(MCIPos.Reconstructed{1,8});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIPos,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIPos.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
BadExecution(1,1) = 1;
BadExecution(4,1) = 1;

%Save the information in the Reconstructed Structure
MCIPos.Reconstructed{1,8}.BadExecution = BadExecution;

%% Participant 9
trialsize = height(MCIPos.Reconstructed{1,9});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIPos,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIPos.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
BadExecution(10,1) = 1;
BadExecution(24,1) = 1;
BadExecution(25,1) = 1;

%Save the information in the Reconstructed Structure
MCIPos.Reconstructed{1,9}.BadExecution = BadExecution;
%% Participant 10
trialsize = height(MCIPos.Reconstructed{1,10});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIPos,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIPos.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
BadExecution(2,1) = 1;
BadExecution(12,1) = 1;
BadExecution(16,1) = 1;

%Save the information in the Reconstructed Structure
MCIPos.Reconstructed{1,10}.BadExecution = BadExecution;
%% Participant 11
trialsize = height(MCIPos.Reconstructed{1,11});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIPos,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIPos.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
BadExecution(2,1) = 1;
BadExecution(3,1) = 1;
BadExecution(4,1) = 1;
BadExecution(5,1) = 1;
BadExecution(8,1) = 1;
BadExecution(9,1) = 1;
BadExecution(15,1) = 1;
BadExecution(16,1) = 1;
BadExecution(17,1) = 1;

%Save the information in the Reconstructed Structure
MCIPos.Reconstructed{1,11}.BadExecution = BadExecution;