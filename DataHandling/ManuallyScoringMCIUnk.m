%Run Calculate TrackingPath first

%% Participant 1
trialsize = height(Unknown.Reconstructed{1,1});
% Prepare structure
BadExecution = zeros(trialsize,1);

% THIS PARTICIPANT MUST BE EXCLUDED. IT DID NOT COMPLETE THE PATH BUT HE TRACED BACK THE TRIANGLE 
BadExecution(1,1) = 1;
BadExecution(2,1) = 1;
BadExecution(3,1) = 1;
BadExecution(4,1) = 1;
BadExecution(5,1) = 1;
BadExecution(6,1) = 1;
BadExecution(7,1) = 1;
BadExecution(8,1) = 1;
BadExecution(9,1) = 1;
BadExecution(10,1) = 1;
BadExecution(11,1) = 1;
BadExecution(12,1) = 1;
BadExecution(13,1) = 1;
BadExecution(14,1) = 1;
BadExecution(15,1) = 1;
BadExecution(16,1) = 1;
BadExecution(17,1) = 1;
BadExecution(18,1) = 1;
BadExecution(19,1) = 1;
BadExecution(20,1) = 1;
BadExecution(21,1) = 1;
BadExecution(22,1) = 1;
BadExecution(23,1) = 1;
BadExecution(24,1) = 1;
BadExecution(25,1) = 1;
BadExecution(26,1) = 1;
BadExecution(27,1) = 1;

%Save the information in the Reconstructed Structure
Unknown.Reconstructed{1,1}.BadExecution = BadExecution;

%% Participant 2
trialsize = height(Unknown.Reconstructed{1,2});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(Unknown,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model Unknown.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection

%Save the information in the Reconstructed Structure
Unknown.Reconstructed{1,2}.BadExecution = BadExecution;

%% Participant 3
trialsize = height(Unknown.Reconstructed{1,3});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(Unknown,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model Unknown.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection

%Save the information in the Reconstructed Structure
Unknown.Reconstructed{1,3}.BadExecution = BadExecution;

%% Participant 4
trialsize = height(Unknown.Reconstructed{1,4});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(Unknown,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model Unknown.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection

%Save the information in the Reconstructed Structure
Unknown.Reconstructed{1,4}.BadExecution = BadExecution;

%% Participant 5
trialsize = height(Unknown.Reconstructed{1,5});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(Unknown,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model Unknown.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
BadExecution(3,1) = 1;
BadExecution(6,1) = 1;
BadExecution(14,1) = 1;
BadExecution(16,1) = 1;
BadExecution(26,1) = 1;
%Save the information in the Reconstructed Structure
Unknown.Reconstructed{1,5}.BadExecution = BadExecution;

%% Participant 6
trialsize = height(Unknown.Reconstructed{1,6});
% Prepare structure
BadExecution = zeros(trialsize,1);

% THIS PARTICIPANT MUST BE EXCLUDED. IT DID NOT FOLLOW THE CONES WHEN
% PROMPTING THEM. DID NOT ATTEMPT TO GIVE AN ANSWER
BadExecution(1,1) = 1;
BadExecution(2,1) = 1;
BadExecution(5,1) = 1;
BadExecution(6,1) = 1;
BadExecution(7,1) = 1;
BadExecution(8,1) = 1;
BadExecution(9,1) = 1;
BadExecution(12,1) = 1;
BadExecution(15,1) = 1;
BadExecution(16,1) = 1;
BadExecution(17,1) = 1;
BadExecution(19,1) = 1;
BadExecution(23,1) = 1;
BadExecution(26,1) = 1;
BadExecution(27,1) = 1;
%Save the information in the Reconstructed Structure
Unknown.Reconstructed{1,6}.BadExecution = BadExecution;

%% Participant 7
trialsize = height(Unknown.Reconstructed{1,7});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Good execution but tracking gives different result
Unknown.Reconstructed{1,7}.RealReturnAngle(2) = Unknown.Reconstructed{1,7}.InferredReturnAngle(2);
Unknown.Reconstructed{1,7}.RealReturnAngle(23) = Unknown.Reconstructed{1,7}.InferredReturnAngle(23);

%Save the information in the Reconstructed Structure
Unknown.Reconstructed{1,7}.BadExecution = BadExecution;

%% Participant 8
trialsize = height(Unknown.Reconstructed{1,8});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(Unknown,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model Unknown.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection

%Save the information in the Reconstructed Structure
Unknown.Reconstructed{1,8}.BadExecution = BadExecution;

%% Participant 9
trialsize = height(Unknown.Reconstructed{1,9});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(Unknown,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model Unknown.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection

%Save the information in the Reconstructed Structure
Unknown.Reconstructed{1,9}.BadExecution = BadExecution;
%% Participant 10
trialsize = height(Unknown.Reconstructed{1,10});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(Unknown,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model Unknown.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection

%Save the information in the Reconstructed Structure
Unknown.Reconstructed{1,10}.BadExecution = BadExecution;
%% Participant 11
trialsize = height(Unknown.Reconstructed{1,11});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(Unknown,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model Unknown.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection

%Save the information in the Reconstructed Structure
Unknown.Reconstructed{1,11}.BadExecution = BadExecution;