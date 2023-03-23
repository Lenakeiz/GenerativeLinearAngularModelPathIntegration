%% Manually scoring the participants trajectory
% In order to the model to fit the data the participant had to demonstrate
% a clear intention of completing the triangle.
% This did not happen all of the time and therefore we had to exclude some
% data. Also we noticed by looking at the tracking data that participants
% turned over the long angle instead of the short one, but calculation
% based on the solely position could not capture this rotation. In order to
% account for real angular production error we corrected the calculated
% angle based on the fact the participant took the long or the short angle
% to complete their estimation.

% Realtime Tracking paths are visualized using
% "VisualizeRealtimeTrackingData" function. Must be run after
% CalculateTrackingPath function.

% Bad trials are marked as follows: 
% Participant did not show a clear intention of moving from cone 3 (we took
% 0.5m as a threshold).
% Participant retraced previous cone to find the position of cone 1.

%% Participant 1
trialsize = height(MCIUnk.Reconstructed{1,1});
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
MCIUnk.Reconstructed{1,1}.BadExecution = BadExecution;

%% Participant 2
trialsize = height(MCIUnk.Reconstructed{1,2});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIUnk,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIUnk.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection

%Save the information in the Reconstructed Structure
MCIUnk.Reconstructed{1,2}.BadExecution = BadExecution;

%% Participant 3
trialsize = height(MCIUnk.Reconstructed{1,3});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIUnk,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIUnk.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection

%Save the information in the Reconstructed Structure
MCIUnk.Reconstructed{1,3}.BadExecution = BadExecution;

%% Participant 4
trialsize = height(MCIUnk.Reconstructed{1,4});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIUnk,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIUnk.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection

%Save the information in the Reconstructed Structure
MCIUnk.Reconstructed{1,4}.BadExecution = BadExecution;

%% Participant 5
trialsize = height(MCIUnk.Reconstructed{1,5});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIUnk,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIUnk.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
BadExecution(3,1) = 1;
BadExecution(6,1) = 1;
BadExecution(14,1) = 1;
BadExecution(16,1) = 1;
BadExecution(26,1) = 1;
%Save the information in the Reconstructed Structure
MCIUnk.Reconstructed{1,5}.BadExecution = BadExecution;

%% Participant 6
trialsize = height(MCIUnk.Reconstructed{1,6});
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
MCIUnk.Reconstructed{1,6}.BadExecution = BadExecution;

%% Participant 7
trialsize = height(MCIUnk.Reconstructed{1,7});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Good execution but tracking gives different result
MCIUnk.Reconstructed{1,7}.RealReturnAngle(2) = MCIUnk.Reconstructed{1,7}.InferredReturnAngle(2);
MCIUnk.Reconstructed{1,7}.RealReturnAngle(23) = MCIUnk.Reconstructed{1,7}.InferredReturnAngle(23);

%Save the information in the Reconstructed Structure
MCIUnk.Reconstructed{1,7}.BadExecution = BadExecution;

%% Participant 8
trialsize = height(MCIUnk.Reconstructed{1,8});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIUnk,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIUnk.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
BadExecution(8,1) = 1;
BadExecution(25,1) = 1;
BadExecution(26,1) = 1;

%Save the information in the Reconstructed Structure
MCIUnk.Reconstructed{1,8}.BadExecution = BadExecution;

%% Participant 9
trialsize = height(MCIUnk.Reconstructed{1,9});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIUnk,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIUnk.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection

%Save the information in the Reconstructed Structure
MCIUnk.Reconstructed{1,9}.BadExecution = BadExecution;
%% Participant 10
trialsize = height(MCIUnk.Reconstructed{1,10});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIUnk,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIUnk.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection

%Save the information in the Reconstructed Structure
MCIUnk.Reconstructed{1,10}.BadExecution = BadExecution;
%% Participant 11
trialsize = height(MCIUnk.Reconstructed{1,11});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIUnk,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIUnk.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
MCIUnk.Reconstructed{1,11}.RealReturnAngle(5) = MCIUnk.Reconstructed{1,11}.InferredReturnAngle(5);

%Save the information in the Reconstructed Structure
MCIUnk.Reconstructed{1,11}.BadExecution = BadExecution;
%% Participant 12
trialsize = height(MCIUnk.Reconstructed{1,12});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIUnk,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIUnk.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection

%Save the information in the Reconstructed Structure
MCIUnk.Reconstructed{1,12}.BadExecution = BadExecution;
%% Participant 13
trialsize = height(MCIUnk.Reconstructed{1,13});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIUnk,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIUnk.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection

%Save the information in the Reconstructed Structure
MCIUnk.Reconstructed{1,13}.BadExecution = BadExecution;
%% Participant 14
trialsize = height(MCIUnk.Reconstructed{1,14});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIUnk,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIUnk.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
MCIUnk.Reconstructed{1,14}.RealReturnAngle(4) = MCIUnk.Reconstructed{1,14}.RealReturnAngle(4) - 360;
MCIUnk.Reconstructed{1,14}.RealReturnAngle(26) = MCIUnk.Reconstructed{1,14}.RealReturnAngle(26) + 360;

%Save the information in the Reconstructed Structure
MCIUnk.Reconstructed{1,14}.BadExecution = BadExecution;
%% Participant 15
trialsize = height(MCIUnk.Reconstructed{1,15});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIUnk,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIUnk.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
BadExecution(4,1) = 1;
MCIUnk.Reconstructed{1,15}.RealReturnAngle(14) = MCIUnk.Reconstructed{1,15}.RealReturnAngle(14) - 360;

%Save the information in the Reconstructed Structure
MCIUnk.Reconstructed{1,15}.BadExecution = BadExecution;
%% Participant 16
%% THIS PARTICIPANT HAS THE OLD PATHS AS WELL! Check if to EXCLUDE (see trial 1 for example
trialsize = height(MCIUnk.Reconstructed{1,16});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIUnk,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIUnk.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
MCIUnk.Reconstructed{1,16}.RealReturnAngle(14) = MCIUnk.Reconstructed{1,16}.InferredReturnAngle(14);

%Save the information in the Reconstructed Structure
MCIUnk.Reconstructed{1,16}.BadExecution = BadExecution;
%% Participant 17
trialsize = height(MCIUnk.Reconstructed{1,17});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIUnk,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIUnk.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
BadExecution(2,1) = 1;
BadExecution(9,1) = 1;

%Save the information in the Reconstructed Structure
MCIUnk.Reconstructed{1,17}.BadExecution = BadExecution;
%% Participant 18
trialsize = height(MCIUnk.Reconstructed{1,18});
% Prepare structure
BadExecution = zeros(trialsize,1);

% Check each single trial by calling VisualizeRealtimeTrackingData(MCIUnk,1,1,0.0001,'cutconethree',true)
% Check the following quantity as well to see what s going to be input on
% the model MCIUnk.Reconstructed{1,1}.RealReturnAngle
% Write here all the bad trials after visual inspection
MCIUnk.Reconstructed{1,18}.RealReturnAngle(3)  = MCIUnk.Reconstructed{1,18}.RealReturnAngle(3) + 360;
BadExecution(4,1) = 1;
MCIUnk.Reconstructed{1,18}.RealReturnAngle(7)  = MCIUnk.Reconstructed{1,18}.RealReturnAngle(7) + 360;
MCIUnk.Reconstructed{1,18}.RealReturnAngle(8)  = MCIUnk.Reconstructed{1,18}.RealReturnAngle(8) + 360;
MCIUnk.Reconstructed{1,18}.RealReturnAngle(25) = MCIUnk.Reconstructed{1,18}.RealReturnAngle(25) - 360;
MCIUnk.Reconstructed{1,18}.RealReturnAngle(27) = MCIUnk.Reconstructed{1,18}.RealReturnAngle(27) + 360;

%Save the information in the Reconstructed Structure
MCIUnk.Reconstructed{1,18}.BadExecution = BadExecution;

%%
clear trialsize BadExecution