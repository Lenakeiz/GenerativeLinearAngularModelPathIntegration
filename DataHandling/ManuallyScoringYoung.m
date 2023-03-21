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
% "VisualizeRealtimeTrackingData" function.

% Bad trials are marked as follows: 
% Participant did not show a clear intention of moving from cone 3 (we took
% 0.5m as a threshold).
% Participant retraced previous cone to find the position of cone 1.

%% Participant 1
trialsize = height(YoungControls.Reconstructed{1,1});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,1}.BadExecution = BadExecution;

%% Participant 2
trialsize = height(YoungControls.Reconstructed{1,2});
% Prepare structure
BadExecution = zeros(trialsize,1);

YoungControls.Reconstructed{1,2}.RealReturnAngle(15) = YoungControls.Reconstructed{1,2}.RealReturnAngle(15) + 360;
YoungControls.Reconstructed{1,2}.RealReturnAngle(21) = YoungControls.Reconstructed{1,2}.RealReturnAngle(21) + 360;
YoungControls.Reconstructed{1,2}.RealReturnAngle(25) = YoungControls.Reconstructed{1,2}.RealReturnAngle(25) + 360;

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,2}.BadExecution = BadExecution;

%% Participant 3
trialsize = height(YoungControls.Reconstructed{1,3});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(11,1) = 1;

YoungControls.Reconstructed{1,3}.RealReturnAngle(14) = YoungControls.Reconstructed{1,3}.RealReturnAngle(14) + 360;

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,3}.BadExecution = BadExecution;

%% Participant 4
trialsize = height(YoungControls.Reconstructed{1,4});
% Prepare structure
BadExecution = zeros(trialsize,1);

YoungControls.Reconstructed{1,4}.RealReturnAngle(24) = YoungControls.Reconstructed{1,4}.RealReturnAngle(24) - 360;
%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,4}.BadExecution = BadExecution;

%% Participant 5
trialsize = height(YoungControls.Reconstructed{1,5});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,5}.BadExecution = BadExecution;

%% Participant 6
trialsize = height(YoungControls.Reconstructed{1,6});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,6}.BadExecution = BadExecution;

%% Participant 7
trialsize = height(YoungControls.Reconstructed{1,7});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,7}.BadExecution = BadExecution;

%% Participant 8
trialsize = height(YoungControls.Reconstructed{1,8});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,8}.BadExecution = BadExecution;

%% Participant 9
trialsize = height(YoungControls.Reconstructed{1,9});
% Prepare structure
BadExecution = zeros(trialsize,1);

YoungControls.Reconstructed{1,9}.RealReturnAngle(11) = YoungControls.Reconstructed{1,9}.RealReturnAngle(11) + 360;
YoungControls.Reconstructed{1,9}.RealReturnAngle(17) = YoungControls.Reconstructed{1,9}.RealReturnAngle(17) + 360;
YoungControls.Reconstructed{1,9}.RealReturnAngle(24) = YoungControls.Reconstructed{1,9}.RealReturnAngle(24) + 360;
YoungControls.Reconstructed{1,9}.RealReturnAngle(34) = YoungControls.Reconstructed{1,9}.RealReturnAngle(34) + 360;

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,9}.BadExecution = BadExecution;

%% Participant 10
trialsize = height(YoungControls.Reconstructed{1,10});
% Prepare structure
BadExecution = zeros(trialsize,1);

YoungControls.Reconstructed{1,10}.RealReturnAngle(5) = YoungControls.Reconstructed{1,10}.RealReturnAngle(5) + 360;
YoungControls.Reconstructed{1,10}.RealReturnAngle(7) = YoungControls.Reconstructed{1,10}.RealReturnAngle(7) + 360;

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,10}.BadExecution = BadExecution;

%% Participant 11
trialsize = height(YoungControls.Reconstructed{1,11});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,11}.BadExecution = BadExecution;

%% Participant 12
trialsize = height(YoungControls.Reconstructed{1,12});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,12}.BadExecution = BadExecution;

%% Participant 13
trialsize = height(YoungControls.Reconstructed{1,13});
% Prepare structure
BadExecution = zeros(trialsize,1);

YoungControls.Reconstructed{1,13}.RealReturnAngle(33) = YoungControls.Reconstructed{1,13}.RealReturnAngle(33) + 360;
%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,13}.BadExecution = BadExecution;

%% Participant 14
trialsize = height(YoungControls.Reconstructed{1,14});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,14}.BadExecution = BadExecution;

%% Participant 15
trialsize = height(YoungControls.Reconstructed{1,15});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,15}.BadExecution = BadExecution;

%% Participant 16
trialsize = height(YoungControls.Reconstructed{1,16});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,16}.BadExecution = BadExecution;

%% Participant 17
trialsize = height(YoungControls.Reconstructed{1,17});
% Prepare structure
BadExecution = zeros(trialsize,1);

YoungControls.Reconstructed{1,17}.RealReturnAngle(23) = YoungControls.Reconstructed{1,17}.RealReturnAngle(23) + 360;
%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,17}.BadExecution = BadExecution;

%% Participant 18
trialsize = height(YoungControls.Reconstructed{1,18});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,18}.BadExecution = BadExecution;

%% Participant 19
trialsize = height(YoungControls.Reconstructed{1,19});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(7,1) = 1;

YoungControls.Reconstructed{1,19}.RealReturnAngle(3) = YoungControls.Reconstructed{1,19}.RealReturnAngle(3) + 360;
%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,19}.BadExecution = BadExecution;

%% Participant 20
trialsize = height(YoungControls.Reconstructed{1,20});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,20}.BadExecution = BadExecution;

%% Participant 21
trialsize = height(YoungControls.Reconstructed{1,21});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,21}.BadExecution = BadExecution;

%% Participant 22
trialsize = height(YoungControls.Reconstructed{1,22});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,22}.BadExecution = BadExecution;

%% Participant 23
trialsize = height(YoungControls.Reconstructed{1,23});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,23}.BadExecution = BadExecution;

%% Participant 24
trialsize = height(YoungControls.Reconstructed{1,24});
% Prepare structure
BadExecution = zeros(trialsize,1);

YoungControls.Reconstructed{1,24}.RealReturnAngle(1) = YoungControls.Reconstructed{1,24}.RealReturnAngle(1) + 360;
YoungControls.Reconstructed{1,24}.RealReturnAngle(4) = YoungControls.Reconstructed{1,24}.RealReturnAngle(4) + 360;
%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,24}.BadExecution = BadExecution;

%% Participant 25
trialsize = height(YoungControls.Reconstructed{1,25});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,25}.BadExecution = BadExecution;

%% Participant 26
trialsize = height(YoungControls.Reconstructed{1,26});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,26}.BadExecution = BadExecution;

%% Participant 27
trialsize = height(YoungControls.Reconstructed{1,27});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(7,1) = 1;

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,27}.BadExecution = BadExecution;

%% Participant 28
trialsize = height(YoungControls.Reconstructed{1,28});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,28}.BadExecution = BadExecution;

%% Participant 29
trialsize = height(YoungControls.Reconstructed{1,29});
% Prepare structure
BadExecution = zeros(trialsize,1);

YoungControls.Reconstructed{1,29}.RealReturnAngle(31) = YoungControls.Reconstructed{1,29}.RealReturnAngle(31) + 360;
%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,29}.BadExecution = BadExecution;

%% Participant 30
trialsize = height(YoungControls.Reconstructed{1,30});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(17,1) = 1;
BadExecution(23,1) = 1;

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,30}.BadExecution = BadExecution;

%% Participant 31
trialsize = height(YoungControls.Reconstructed{1,31});
% Prepare structure
BadExecution = zeros(trialsize,1);

YoungControls.Reconstructed{1,31}.RealReturnAngle(12) = YoungControls.Reconstructed{1,31}.RealReturnAngle(12) + 360;
YoungControls.Reconstructed{1,31}.RealReturnAngle(13) = YoungControls.Reconstructed{1,31}.RealReturnAngle(13) + 360;
YoungControls.Reconstructed{1,31}.RealReturnAngle(26) = YoungControls.Reconstructed{1,31}.RealReturnAngle(26) + 360;
YoungControls.Reconstructed{1,31}.RealReturnAngle(35) = YoungControls.Reconstructed{1,31}.RealReturnAngle(35) + 360;

%Save the information in the Reconstructed Structure
YoungControls.Reconstructed{1,31}.BadExecution = BadExecution;

%%
clear trialsize BadExecution