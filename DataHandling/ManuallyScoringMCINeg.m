%Run Calculate TrackingPath first

%% Participant 1
trialsize = height(MCINeg.Reconstructed{1,1});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCINeg,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIPos.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection

%Save the information in the Reconstructed Structure
MCINeg.Reconstructed{1,1}.BadExecution = BadExecution;

%% Participant 2
trialsize = height(MCINeg.Reconstructed{1,2});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(1,1) = 1;

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIPos,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIPos.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection

%Save the information in the Reconstructed Structure
MCINeg.Reconstructed{1,2}.BadExecution = BadExecution;

%% Participant 3
trialsize = height(MCINeg.Reconstructed{1,3});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(1,1) = 1;
BadExecution(2,1) = 1;

%Save the information in the Reconstructed Structure
MCINeg.Reconstructed{1,3}.BadExecution = BadExecution;

%% Participant 4
trialsize = height(MCINeg.Reconstructed{1,4});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(8,1) = 1;
BadExecution(10,1) = 1;
BadExecution(12,1) = 1;

MCINeg.Reconstructed{1,4}.RealReturnAngle(7) = MCINeg.Reconstructed{1,4}.RealReturnAngle(7) + 360;

%Save the information in the Reconstructed Structure
MCINeg.Reconstructed{1,4}.BadExecution = BadExecution;

%% Participant 5
trialsize = height(MCINeg.Reconstructed{1,5});
% Prepare structure
BadExecution = zeros(trialsize,1);

MCINeg.Reconstructed{1,5}.RealReturnAngle(12) = MCINeg.Reconstructed{1,5}.RealReturnAngle(12) + 360;

%Save the information in the Reconstructed Structure
MCINeg.Reconstructed{1,5}.BadExecution = BadExecution;

%% Participant 6
trialsize = height(MCINeg.Reconstructed{1,6});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(2,1) = 1;
BadExecution(15,1) = 1;

%Save the information in the Reconstructed Structure
MCINeg.Reconstructed{1,6}.BadExecution = BadExecution;

%% Participant 7
trialsize = height(MCINeg.Reconstructed{1,7});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(1,1) = 1;

%Save the information in the Reconstructed Structure
MCINeg.Reconstructed{1,7}.BadExecution = BadExecution;

%% Participant 8
trialsize = height(MCINeg.Reconstructed{1,8});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(1,1) = 1;

MCINeg.Reconstructed{1,8}.RealReturnAngle(1) = MCINeg.Reconstructed{1,8}.RealReturnAngle(1) - 360;

%Save the information in the Reconstructed Structure
MCINeg.Reconstructed{1,8}.BadExecution = BadExecution;

%% Participant 9
trialsize = height(MCINeg.Reconstructed{1,9});
% Prepare structure
BadExecution = zeros(trialsize,1);

MCINeg.Reconstructed{1,9}.RealReturnAngle(6) = MCINeg.Reconstructed{1,9}.RealReturnAngle(6) + 360;
MCINeg.Reconstructed{1,9}.RealReturnAngle(17) = MCINeg.Reconstructed{1,9}.RealReturnAngle(17) + 360;

%Save the information in the Reconstructed Structure
MCINeg.Reconstructed{1,9}.BadExecution = BadExecution;

%% Participant 10
trialsize = height(MCINeg.Reconstructed{1,10});
% Prepare structure
BadExecution = zeros(trialsize,1);


%Save the information in the Reconstructed Structure
MCINeg.Reconstructed{1,10}.BadExecution = BadExecution;