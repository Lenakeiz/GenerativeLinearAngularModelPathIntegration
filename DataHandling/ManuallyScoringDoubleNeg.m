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
trialsize = height(DoubleNeg.Reconstructed{1,1});
% Prepare structure
BadExecution = zeros(trialsize,1);

DoubleNeg.Reconstructed{1,1}.RealReturnAngle(26) = DoubleNeg.Reconstructed{1,1}.RealReturnAngle(26) + 360;
DoubleNeg.Reconstructed{1,1}.RealReturnAngle(27) = DoubleNeg.Reconstructed{1,1}.RealReturnAngle(27) + 360;
DoubleNeg.Reconstructed{1,1}.RealReturnAngle(33) = DoubleNeg.Reconstructed{1,1}.RealReturnAngle(33) + 360;

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,1}.BadExecution = BadExecution;

%% Participant 2
trialsize = height(DoubleNeg.Reconstructed{1,2});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,2}.BadExecution = BadExecution;

%% Participant 3
trialsize = height(DoubleNeg.Reconstructed{1,3});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(11,1) = 1;

DoubleNeg.Reconstructed{1,3}.RealReturnAngle(4) = DoubleNeg.Reconstructed{1,3}.RealReturnAngle(4) + 360;
DoubleNeg.Reconstructed{1,3}.RealReturnAngle(5) = DoubleNeg.Reconstructed{1,3}.RealReturnAngle(5) + 360;
DoubleNeg.Reconstructed{1,3}.RealReturnAngle(24) = DoubleNeg.Reconstructed{1,3}.RealReturnAngle(24) + 360;

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,3}.BadExecution = BadExecution;

%% Participant 4
trialsize = height(DoubleNeg.Reconstructed{1,4});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,4}.BadExecution = BadExecution;

%% Participant 5
trialsize = height(DoubleNeg.Reconstructed{1,5});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,5}.BadExecution = BadExecution;

%% Participant 6
trialsize = height(DoubleNeg.Reconstructed{1,6});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,6}.BadExecution = BadExecution;

%% Participant 7
trialsize = height(DoubleNeg.Reconstructed{1,7});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,7}.BadExecution = BadExecution;

%% Participant 8
trialsize = height(DoubleNeg.Reconstructed{1,8});
% Prepare structure
BadExecution = zeros(trialsize,1);

DoubleNeg.Reconstructed{1,8}.RealReturnAngle(4) = DoubleNeg.Reconstructed{1,8}.RealReturnAngle(4) + 360;
DoubleNeg.Reconstructed{1,8}.RealReturnAngle(7) = DoubleNeg.Reconstructed{1,8}.RealReturnAngle(7) + 360;
DoubleNeg.Reconstructed{1,8}.RealReturnAngle(8) = DoubleNeg.Reconstructed{1,8}.RealReturnAngle(8) + 360;
DoubleNeg.Reconstructed{1,8}.RealReturnAngle(11) = DoubleNeg.Reconstructed{1,8}.RealReturnAngle(11) + 360;
DoubleNeg.Reconstructed{1,8}.RealReturnAngle(12) = DoubleNeg.Reconstructed{1,8}.RealReturnAngle(12) + 360;
DoubleNeg.Reconstructed{1,8}.RealReturnAngle(14) = DoubleNeg.Reconstructed{1,8}.RealReturnAngle(14) - 360;
DoubleNeg.Reconstructed{1,8}.RealReturnAngle(17) = DoubleNeg.Reconstructed{1,8}.RealReturnAngle(17) + 360;
DoubleNeg.Reconstructed{1,8}.RealReturnAngle(24) = DoubleNeg.Reconstructed{1,8}.RealReturnAngle(24) + 360;

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,8}.BadExecution = BadExecution;

%% Participant 9
trialsize = height(DoubleNeg.Reconstructed{1,9});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(32,1) = 1;

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,9}.BadExecution = BadExecution;


%% Participant 10
trialsize = height(DoubleNeg.Reconstructed{1,10});
% Prepare structure
BadExecution = zeros(trialsize,1);

DoubleNeg.Reconstructed{1,10}.RealReturnAngle(4) = DoubleNeg.Reconstructed{1,10}.RealReturnAngle(4) - 360;

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,10}.BadExecution = BadExecution;

%% Participant 11
trialsize = height(DoubleNeg.Reconstructed{1,11});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,11}.BadExecution = BadExecution;

%% Participant 12
trialsize = height(DoubleNeg.Reconstructed{1,12});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,12}.BadExecution = BadExecution;

%% Participant 13
trialsize = height(DoubleNeg.Reconstructed{1,13});
% Prepare structure
BadExecution = zeros(trialsize,1);

DoubleNeg.Reconstructed{1,13}.RealReturnAngle(2) = DoubleNeg.Reconstructed{1,13}.RealReturnAngle(2) - 360;


%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,13}.BadExecution = BadExecution;

%% Participant 14
trialsize = height(DoubleNeg.Reconstructed{1,14});
% Prepare structure
BadExecution = zeros(trialsize,1);

DoubleNeg.Reconstructed{1,14}.RealReturnAngle(27) = DoubleNeg.Reconstructed{1,14}.RealReturnAngle(27) + 360;


%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,14}.BadExecution = BadExecution;

%% Participant 15
trialsize = height(DoubleNeg.Reconstructed{1,15});
% Prepare structure
BadExecution = zeros(trialsize,1);

DoubleNeg.Reconstructed{1,15}.RealReturnAngle(19) = DoubleNeg.Reconstructed{1,15}.RealReturnAngle(19) + 360;
DoubleNeg.Reconstructed{1,15}.RealReturnAngle(24) = DoubleNeg.Reconstructed{1,15}.RealReturnAngle(24) + 360;
DoubleNeg.Reconstructed{1,15}.RealReturnAngle(29) = DoubleNeg.Reconstructed{1,15}.RealReturnAngle(29) + 360;

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,15}.BadExecution = BadExecution;


%% Participant 16
trialsize = height(DoubleNeg.Reconstructed{1,16});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,16}.BadExecution = BadExecution;

%% Participant 17
trialsize = height(DoubleNeg.Reconstructed{1,17});
% Prepare structure
BadExecution = zeros(trialsize,1);

DoubleNeg.Reconstructed{1,17}.RealReturnAngle(22) = DoubleNeg.Reconstructed{1,17}.RealReturnAngle(22) + 360;
DoubleNeg.Reconstructed{1,17}.RealReturnAngle(28) = DoubleNeg.Reconstructed{1,17}.RealReturnAngle(28) + 360;

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,17}.BadExecution = BadExecution;

%% Participant 18
trialsize = height(DoubleNeg.Reconstructed{1,18});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(7,1) = 1;

DoubleNeg.Reconstructed{1,18}.RealReturnAngle(24) = DoubleNeg.Reconstructed{1,18}.RealReturnAngle(24) + 360;

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,18}.BadExecution = BadExecution;

%% Participant 19
trialsize = height(DoubleNeg.Reconstructed{1,19});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,19}.BadExecution = BadExecution;

%% Participant 20
trialsize = height(DoubleNeg.Reconstructed{1,20});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,20}.BadExecution = BadExecution;

%% Participant 21
trialsize = height(DoubleNeg.Reconstructed{1,21});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,21}.BadExecution = BadExecution;

%% Participant 22
trialsize = height(DoubleNeg.Reconstructed{1,22});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,22}.BadExecution = BadExecution;

%% Participant 23
trialsize = height(DoubleNeg.Reconstructed{1,23});
% Prepare structure
BadExecution = zeros(trialsize,1);

DoubleNeg.Reconstructed{1,23}.RealReturnAngle(5) = DoubleNeg.Reconstructed{1,23}.RealReturnAngle(5) + 360;
DoubleNeg.Reconstructed{1,23}.RealReturnAngle(6) = DoubleNeg.Reconstructed{1,23}.RealReturnAngle(6) + 360;


%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,23}.BadExecution = BadExecution;

%% Participant 24

trialsize = height(DoubleNeg.Reconstructed{1,24});
% Prepare structure
BadExecution = zeros(trialsize,1);

DoubleNeg.Reconstructed{1,24}.RealReturnAngle(12) = DoubleNeg.Reconstructed{1,24}.RealReturnAngle(12) + 360;
DoubleNeg.Reconstructed{1,24}.RealReturnAngle(25) = DoubleNeg.Reconstructed{1,24}.RealReturnAngle(25) + 360;
DoubleNeg.Reconstructed{1,24}.RealReturnAngle(30) = DoubleNeg.Reconstructed{1,24}.RealReturnAngle(30) + 360;


%Save the information in the Reconstructed Structure
DoubleNeg.Reconstructed{1,24}.BadExecution = BadExecution;

%%
clear trialsize BadExecution