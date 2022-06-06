%Run Calculate TrackingPath first


% the goal of this script is to manually identify ' bad trials'
% in participants from the healthy old group. The judgement is subjective
% and is based on observing the participant's trajectory - and judging
% whether it looks like they know what they are doing. 

% if trial is judged as a bad trial, then mark it as follows:
% BadExecution(trial number, 1) = 1;
% otherwise BadExecution(trial number, 1) = 0 as it was originally
% initialised to 0. 
% note that BadExecution(trial number, 1) has 1 in the second variable, and
% that is always the case because that represents there is only one column
% participant number associated with given BadExecution are assigned at the
% end of each section, when we save the information in the structure
% corresponding to each participant separately 

% note:  
% MCINeg.Reconstructed{1,participant number}.RealReturnAngle(trial number) =
% MCINeg.Reconstructed{1,participant number}.RealReturnAngle(trial number) + 360;

% use this line of code, when you identify that the 'Real Inferred Angle'
% stated in the table in the visualisation has a 'weird' angle 
% (usually a weirdly large negative angle) - this will be the case when 
% we erroneously deducted (-360 degrees) 
% we need to add the 360 degrees back. 
% Note that originally we deducted 360 degrees supposedely in trials, where
% participants turned in the wrong ('clockwise') direction at cone 3 just
% prior to embarking on the return journey. If that really happened - i.e.
% if at cone 3 just before outbound path they truly turned to the clockwise
% direction, then 360 degrees SHOULD be deducted. This would result in the
% first 'body rotation' value being negative, reflecting the clockwise
% direction of turning (when anticlockwise was expected). Since the 'body
% rotation' would be negative, the 'inferred angle' should be positive, 
% these two angles would have opposing signs, and
% altogether the 'real inferred angle' would have 360 deducted

% However, if we visually see that the 'real inferred angle' did have the
% 360 degrees deducted, and we see that the 'body rotation' angle is
% negative - * but not because they initially turned clockwise prior to
% returning, but because they later in the trajectory turned around
% clockwise, which altogether caused the negative value * - then we should
% add the 360 degrees back, as it was a mistake to deduce it to begin with 

% This error happens when the 'body rotation' and 'inferred angle' were of
% different signs, and hence (-360) degrees were deducted, but if we
% manually look at the trajectory, it should not have been deducted, so we
% add it back 

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

BadExecution(10,1) = 1;


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

BadExecution(25,1) = 1;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,3}.BadExecution = BadExecution;

%% Participant 4
trialsize = height(HealthyControls.Reconstructed{1,4});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(7,1) = 1;

HealthyControls.Reconstructed{1,4}.RealReturnAngle(5) = HealthyControls.Reconstructed{1,4}.RealReturnAngle(5) + 360;
HealthyControls.Reconstructed{1,4}.RealReturnAngle(11) = HealthyControls.Reconstructed{1,4}.RealReturnAngle(11) + 360;
HealthyControls.Reconstructed{1,4}.RealReturnAngle(15) = HealthyControls.Reconstructed{1,4}.RealReturnAngle(15) + 360;
HealthyControls.Reconstructed{1,4}.RealReturnAngle(18) = HealthyControls.Reconstructed{1,4}.RealReturnAngle(18) + 360;

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
HealthyControls.Reconstructed{1,9}.RealReturnAngle(18) = HealthyControls.Reconstructed{1,9}.RealReturnAngle(18) + 360;


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
BadExecution(12,1) = 1;

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

HealthyControls.Reconstructed{1,18}.RealReturnAngle(25) = HealthyControls.Reconstructed{1,18}.RealReturnAngle(25) + 360;

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
BadExecution(16,1) = 1;
BadExecution(21,1) = 1;
BadExecution(21,1) = 1;

HealthyControls.Reconstructed{1,26}.RealReturnAngle(4) = HealthyControls.Reconstructed{1,26}.RealReturnAngle(4) + 360;
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

HealthyControls.Reconstructed{1,27}.RealReturnAngle(8) = HealthyControls.Reconstructed{1,27}.RealReturnAngle(8) + 360;
% double check trial 9 - walking backwards??
HealthyControls.Reconstructed{1,26}.RealReturnAngle(9) = HealthyControls.Reconstructed{1,26}.RealReturnAngle(9) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,27}.BadExecution = BadExecution;

%% Participant 28
trialsize = height(HealthyControls.Reconstructed{1,28});
% Prepare structure
BadExecution = zeros(trialsize,1);

HealthyControls.Reconstructed{1,28}.RealReturnAngle(21) = HealthyControls.Reconstructed{1,28}.RealReturnAngle(21) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,28}.BadExecution = BadExecution;

%% Participant 29
trialsize = height(HealthyControls.Reconstructed{1,29});
% Prepare structure
BadExecution = zeros(trialsize,1);

%ATTENTION - DOUBLE CHECK TRIAL 14
BadExecution(14,1) = 1;

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
HealthyControls.Reconstructed{1,32}.RealReturnAngle(18) = HealthyControls.Reconstructed{1,32}.RealReturnAngle(18) + 360;
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