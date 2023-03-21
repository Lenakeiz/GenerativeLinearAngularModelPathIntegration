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
trialsize = height(HealthyControls.Reconstructed{1,1});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,1}.BadExecution = BadExecution;

%% Participant 2
trialsize = height(HealthyControls.Reconstructed{1,2});
% Prepare structure
BadExecution = zeros(trialsize,1);

% real inferred angle - if weird, add back +360 degrees. this is when the
% value is suspisciusly large (as -360 were deducted). If a trial is
% already a bad trial, then no need for this action. Note: a trial may be a
% good trial but may still require adding back the 360 degrees   
HealthyControls.Reconstructed{1,2}.RealReturnAngle(20) = HealthyControls.Reconstructed{1,2}.RealReturnAngle(20) + 360;
HealthyControls.Reconstructed{1,2}.RealReturnAngle(23) = HealthyControls.Reconstructed{1,2}.RealReturnAngle(23) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,2}.BadExecution = BadExecution;

%% Participant 3
trialsize = height(HealthyControls.Reconstructed{1,3});
% Prepare structure
BadExecution = zeros(trialsize,1);

%BadExecution(25,1) = 1; %???????????????????

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,3}.BadExecution = BadExecution;

%% Participant 4
trialsize = height(HealthyControls.Reconstructed{1,4});
% Prepare structure
BadExecution = zeros(trialsize,1);

%BadExecution(7,1) = 1;

HealthyControls.Reconstructed{1,4}.RealReturnAngle(5) = HealthyControls.Reconstructed{1,4}.RealReturnAngle(5) + 360;
HealthyControls.Reconstructed{1,4}.RealReturnAngle(11) = HealthyControls.Reconstructed{1,4}.RealReturnAngle(11) + 360;
HealthyControls.Reconstructed{1,4}.RealReturnAngle(15) = HealthyControls.Reconstructed{1,4}.RealReturnAngle(15) + 360;
%HealthyControls.Reconstructed{1,4}.RealReturnAngle(18) = HealthyControls.Reconstructed{1,4}.RealReturnAngle(18) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,4}.BadExecution = BadExecution;

%% Participant 5
trialsize = height(HealthyControls.Reconstructed{1,5});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,5}.BadExecution = BadExecution;

%% Participant 6
trialsize = height(HealthyControls.Reconstructed{1,6});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(1,1) = 1;
% not sure about 11
BadExecution(11,1) = 1;
BadExecution(26,1) = 1;


HealthyControls.Reconstructed{1,6}.RealReturnAngle(25) = HealthyControls.Reconstructed{1,6}.RealReturnAngle(25) + 360;
HealthyControls.Reconstructed{1,6}.RealReturnAngle(27) = HealthyControls.Reconstructed{1,6}.RealReturnAngle(27) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,6}.BadExecution = BadExecution;

%% Participant 7
trialsize = height(HealthyControls.Reconstructed{1,7});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,7}.BadExecution = BadExecution;

%% Participant 8
trialsize = height(HealthyControls.Reconstructed{1,8});
% Prepare structure
BadExecution = zeros(trialsize,1);


%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,8}.BadExecution = BadExecution;

%% Participant 9
trialsize = height(HealthyControls.Reconstructed{1,9});
% Prepare structure
BadExecution = zeros(trialsize,1);

HealthyControls.Reconstructed{1,9}.RealReturnAngle(5) = HealthyControls.Reconstructed{1,9}.RealReturnAngle(5) + 360;
HealthyControls.Reconstructed{1,9}.RealReturnAngle(11) = HealthyControls.Reconstructed{1,9}.RealReturnAngle(11) + 360;
HealthyControls.Reconstructed{1,9}.RealReturnAngle(15) = HealthyControls.Reconstructed{1,9}.RealReturnAngle(15) + 360;
%HealthyControls.Reconstructed{1,9}.RealReturnAngle(18) = HealthyControls.Reconstructed{1,9}.RealReturnAngle(18) + 360;


%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,9}.BadExecution = BadExecution;

%% Participant 10
trialsize = height(HealthyControls.Reconstructed{1,10});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,10}.BadExecution = BadExecution;

%% Participant 11
trialsize = height(HealthyControls.Reconstructed{1,11});
% Prepare structure
BadExecution = zeros(trialsize,1);

% trial 12 - double check why angle is so big 
%BadExecution(12,1) = 1;

HealthyControls.Reconstructed{1,11}.RealReturnAngle(12) = HealthyControls.Reconstructed{1,11}.RealReturnAngle(12) - 360;
HealthyControls.Reconstructed{1,11}.RealReturnAngle(15) = HealthyControls.Reconstructed{1,11}.RealReturnAngle(15) + 360;
HealthyControls.Reconstructed{1,11}.RealReturnAngle(19) = HealthyControls.Reconstructed{1,11}.RealReturnAngle(19) + 360;
HealthyControls.Reconstructed{1,11}.RealReturnAngle(21) = HealthyControls.Reconstructed{1,11}.RealReturnAngle(21) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,11}.BadExecution = BadExecution;

%% Participant 12
trialsize = height(HealthyControls.Reconstructed{1,12});
% Prepare structure
BadExecution = zeros(trialsize,1);

HealthyControls.Reconstructed{1,12}.RealReturnAngle(4) = HealthyControls.Reconstructed{1,12}.RealReturnAngle(4) + 360;
HealthyControls.Reconstructed{1,12}.RealReturnAngle(15) = HealthyControls.Reconstructed{1,12}.RealReturnAngle(15) + 360;
HealthyControls.Reconstructed{1,12}.RealReturnAngle(26) = HealthyControls.Reconstructed{1,12}.RealReturnAngle(26) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,12}.BadExecution = BadExecution;

%% Participant 13
trialsize = height(HealthyControls.Reconstructed{1,13});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,13}.BadExecution = BadExecution;

%% Participant 14
trialsize = height(HealthyControls.Reconstructed{1,14});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,14}.BadExecution = BadExecution;

%% Participant 15
trialsize = height(HealthyControls.Reconstructed{1,15});
% Prepare structure
BadExecution = zeros(trialsize,1);

% ATTENTION: NAN OUT-OF-BOUND TRIAL!!!
BadExecution(22,1) = 1;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,15}.BadExecution = BadExecution;

%% Participant 16
trialsize = height(HealthyControls.Reconstructed{1,16});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,16}.BadExecution = BadExecution;

%% Participant 17
trialsize = height(HealthyControls.Reconstructed{1,17});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,17}.BadExecution = BadExecution;

%% Participant 18
trialsize = height(HealthyControls.Reconstructed{1,18});
% Prepare structure
BadExecution = zeros(trialsize,1);

%HealthyControls.Reconstructed{1,18}.RealReturnAngle(25) = HealthyControls.Reconstructed{1,18}.RealReturnAngle(25) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,18}.BadExecution = BadExecution;

%% Participant 19
trialsize = height(HealthyControls.Reconstructed{1,19});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,19}.BadExecution = BadExecution;

%% Participant 20
trialsize = height(HealthyControls.Reconstructed{1,20});
% Prepare structure
BadExecution = zeros(trialsize,1);

% ATTENTION: NAN OUT-OF-BOUND TRIAL!!!
BadExecution(5,1) = 1;
BadExecution(7,1) = 1;
BadExecution(8,1) = 1;
BadExecution(11,1) = 1;
BadExecution(17,1) = 1;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,20}.BadExecution = BadExecution;

%% Participant 21
trialsize = height(HealthyControls.Reconstructed{1,21});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,21}.BadExecution = BadExecution;

%% Participant 22
trialsize = height(HealthyControls.Reconstructed{1,22});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,22}.BadExecution = BadExecution;

%% Participant 23
trialsize = height(HealthyControls.Reconstructed{1,23});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,23}.BadExecution = BadExecution;

%% Participant 24
trialsize = height(HealthyControls.Reconstructed{1,24});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,24}.BadExecution = BadExecution;

%% Participant 25
trialsize = height(HealthyControls.Reconstructed{1,25});
% Prepare structure
BadExecution = zeros(trialsize,1);

HealthyControls.Reconstructed{1,25}.RealReturnAngle(20) = HealthyControls.Reconstructed{1,25}.RealReturnAngle(20) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,25}.BadExecution = BadExecution;

%% Participant 26
trialsize = height(HealthyControls.Reconstructed{1,26});
% Prepare structure
BadExecution = zeros(trialsize,1);

% double check: is 16 a bad trial (little walking) ad why the real interred
% angle is so large?
%is 21 a bad trial(very little walking)?
%BadExecution(16,1) = 1;
%BadExecution(21,1) = 1;

HealthyControls.Reconstructed{1,26}.RealReturnAngle(4) = HealthyControls.Reconstructed{1,26}.RealReturnAngle(4) + 360;
HealthyControls.Reconstructed{1,26}.RealReturnAngle(16) = HealthyControls.Reconstructed{1,26}.RealReturnAngle(16) - 360;
HealthyControls.Reconstructed{1,26}.RealReturnAngle(23) = HealthyControls.Reconstructed{1,26}.RealReturnAngle(23) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,26}.BadExecution = BadExecution;

%% Participant 27
trialsize = height(HealthyControls.Reconstructed{1,27});
% Prepare structure
BadExecution = zeros(trialsize,1);

% ATTENTION: NAN OUT-OF-BOUND TRIAL!!!
BadExecution(5,1) = 1;
BadExecution(12,1) = 1;
BadExecution(23,1) = 1;
BadExecution(25,1) = 1;

%HealthyControls.Reconstructed{1,27}.RealReturnAngle(8) = HealthyControls.Reconstructed{1,27}.RealReturnAngle(8) + 360;
% double check trial 9 - walking backwards??
%HealthyControls.Reconstructed{1,26}.RealReturnAngle(9) = HealthyControls.Reconstructed{1,26}.RealReturnAngle(9) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,27}.BadExecution = BadExecution;

%% Participant 28
trialsize = height(HealthyControls.Reconstructed{1,28});
% Prepare structure
BadExecution = zeros(trialsize,1);

%HealthyControls.Reconstructed{1,28}.RealReturnAngle(21) = HealthyControls.Reconstructed{1,28}.RealReturnAngle(21) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,28}.BadExecution = BadExecution;

%% Participant 29
trialsize = height(HealthyControls.Reconstructed{1,29});
% Prepare structure
BadExecution = zeros(trialsize,1);

%ATTENTION - DOUBLE CHECK TRIAL 14
%BadExecution(14,1) = 1;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,29}.BadExecution = BadExecution;

%% Participant 30
trialsize = height(HealthyControls.Reconstructed{1,30});
% Prepare structure
BadExecution = zeros(trialsize,1);

HealthyControls.Reconstructed{1,30}.RealReturnAngle(21) = HealthyControls.Reconstructed{1,30}.RealReturnAngle(21) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,30}.BadExecution = BadExecution;

%% Participant 31
trialsize = height(HealthyControls.Reconstructed{1,31});
% Prepare structure
BadExecution = zeros(trialsize,1);

%check trial 6

HealthyControls.Reconstructed{1,31}.RealReturnAngle(1) = HealthyControls.Reconstructed{1,31}.RealReturnAngle(1) + 360;
HealthyControls.Reconstructed{1,31}.RealReturnAngle(6) = HealthyControls.Reconstructed{1,31}.RealReturnAngle(6) - 360;
HealthyControls.Reconstructed{1,31}.RealReturnAngle(20) = HealthyControls.Reconstructed{1,31}.RealReturnAngle(20) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,31}.BadExecution = BadExecution;

%% Participant 32
trialsize = height(HealthyControls.Reconstructed{1,32});
% Prepare structure
BadExecution = zeros(trialsize,1);

% double check trial 10

% double check trial 6
HealthyControls.Reconstructed{1,32}.RealReturnAngle(6) = HealthyControls.Reconstructed{1,32}.RealReturnAngle(6) + 360;
%HealthyControls.Reconstructed{1,32}.RealReturnAngle(18) = HealthyControls.Reconstructed{1,32}.RealReturnAngle(18) + 360;
HealthyControls.Reconstructed{1,32}.RealReturnAngle(19) = HealthyControls.Reconstructed{1,32}.RealReturnAngle(19) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,32}.BadExecution = BadExecution;

%% Participant 33
trialsize = height(HealthyControls.Reconstructed{1,33});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,33}.BadExecution = BadExecution;

%% Participant 34
trialsize = height(HealthyControls.Reconstructed{1,34});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,34}.BadExecution = BadExecution;

%% Participant 35
trialsize = height(HealthyControls.Reconstructed{1,35});
% Prepare structure
BadExecution = zeros(trialsize,1);

HealthyControls.Reconstructed{1,35}.RealReturnAngle(24) = HealthyControls.Reconstructed{1,35}.RealReturnAngle(24) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,35}.BadExecution = BadExecution;

%% Participant 36
trialsize = height(HealthyControls.Reconstructed{1,36});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,36}.BadExecution = BadExecution;

%%
clear trialsize BadExecution