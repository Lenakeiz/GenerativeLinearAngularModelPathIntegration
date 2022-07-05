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
trialsize = height(FamilyHistPos.Reconstructed{1,1});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,1}.RealReturnAngle(15) = FamilyHistPos.Reconstructed{1,1}.RealReturnAngle(15) + 360;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,1}.BadExecution = BadExecution;

%% Participant 2
trialsize = height(FamilyHistPos.Reconstructed{1,2});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(1,1) = 1;
BadExecution(4,1) = 1;
BadExecution(7,1) = 1;
BadExecution(29,1) = 1;

FamilyHistPos.Reconstructed{1,2}.RealReturnAngle(21) = FamilyHistPos.Reconstructed{1,2}.RealReturnAngle(21) + 360;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,2}.BadExecution = BadExecution;

%% Participant 3
trialsize = height(FamilyHistPos.Reconstructed{1,3});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,3}.RealReturnAngle(10) = FamilyHistPos.Reconstructed{1,3}.RealReturnAngle(10) + 360;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,3}.BadExecution = BadExecution;

%% Participant 4
trialsize = height(FamilyHistPos.Reconstructed{1,4});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,4}.BadExecution = BadExecution;

%% Participant 5
trialsize = height(FamilyHistPos.Reconstructed{1,5});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,5}.BadExecution = BadExecution;

%% Participant 6
trialsize = height(FamilyHistPos.Reconstructed{1,6});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,6}.RealReturnAngle(8) = FamilyHistPos.Reconstructed{1,6}.RealReturnAngle(8) + 360;
FamilyHistPos.Reconstructed{1,6}.RealReturnAngle(14) = FamilyHistPos.Reconstructed{1,6}.RealReturnAngle(14) + 360;
FamilyHistPos.Reconstructed{1,6}.RealReturnAngle(19) = FamilyHistPos.Reconstructed{1,6}.RealReturnAngle(19) + 360;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,6}.BadExecution = BadExecution;

%% Participant 7
trialsize = height(FamilyHistPos.Reconstructed{1,7});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,7}.BadExecution = BadExecution;

%% Participant 8
trialsize = height(FamilyHistPos.Reconstructed{1,8});
% Prepare structure
BadExecution = zeros(trialsize,1);


%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,8}.BadExecution = BadExecution;

%% Participant 9
trialsize = height(FamilyHistPos.Reconstructed{1,9});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,9}.BadExecution = BadExecution;

%% Participant 10
trialsize = height(FamilyHistPos.Reconstructed{1,10});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,10}.BadExecution = BadExecution;

%% Participant 11
trialsize = height(FamilyHistPos.Reconstructed{1,11});
% Prepare structure
BadExecution = zeros(trialsize,1);

% trial 12 - double check why angle is so big 
%BadExecution(12,1) = 1;

FamilyHistPos.Reconstructed{1,11}.RealReturnAngle(3) = FamilyHistPos.Reconstructed{1,11}.RealReturnAngle(3) + 360;
FamilyHistPos.Reconstructed{1,11}.RealReturnAngle(28) = FamilyHistPos.Reconstructed{1,11}.RealReturnAngle(28) + 360;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,11}.BadExecution = BadExecution;

%% Participant 12
trialsize = height(FamilyHistPos.Reconstructed{1,12});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(4,1) = 1;
BadExecution(7,1) = 1;
BadExecution(8,1) = 1;
BadExecution(12,1) = 1;
BadExecution(15,1) = 1;
BadExecution(18,1) = 1;
BadExecution(31,1) = 1;

FamilyHistPos.Reconstructed{1,12}.RealReturnAngle(2) = FamilyHistPos.Reconstructed{1,12}.RealReturnAngle(2) + 360;
FamilyHistPos.Reconstructed{1,12}.RealReturnAngle(9) = FamilyHistPos.Reconstructed{1,12}.RealReturnAngle(9) + 360;
FamilyHistPos.Reconstructed{1,12}.RealReturnAngle(11) = FamilyHistPos.Reconstructed{1,12}.RealReturnAngle(11) + 360;
FamilyHistPos.Reconstructed{1,12}.RealReturnAngle(30) = FamilyHistPos.Reconstructed{1,12}.RealReturnAngle(30) + 360;
FamilyHistPos.Reconstructed{1,12}.RealReturnAngle(33) = FamilyHistPos.Reconstructed{1,12}.RealReturnAngle(33) + 360;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,12}.BadExecution = BadExecution;

%% Participant 13
trialsize = height(FamilyHistPos.Reconstructed{1,13});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,13}.BadExecution = BadExecution;

%% Participant 14
trialsize = height(FamilyHistPos.Reconstructed{1,14});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,14}.BadExecution = BadExecution;

%% Participant 15
trialsize = height(FamilyHistPos.Reconstructed{1,15});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,15}.BadExecution = BadExecution;

%% Participant 16
trialsize = height(FamilyHistPos.Reconstructed{1,16});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,16}.BadExecution = BadExecution;

%% Participant 17
trialsize = height(FamilyHistPos.Reconstructed{1,17});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,17}.BadExecution = BadExecution;

%% Participant 18
trialsize = height(FamilyHistPos.Reconstructed{1,18});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,18}.BadExecution = BadExecution;

%% Participant 19
trialsize = height(FamilyHistPos.Reconstructed{1,19});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,19}.RealReturnAngle(6) = FamilyHistPos.Reconstructed{1,19}.RealReturnAngle(6) + 360;
FamilyHistPos.Reconstructed{1,19}.RealReturnAngle(11) = FamilyHistPos.Reconstructed{1,19}.RealReturnAngle(11) + 360;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,19}.BadExecution = BadExecution;

%% Participant 20
trialsize = height(FamilyHistPos.Reconstructed{1,20});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,20}.BadExecution = BadExecution;

%% Participant 21
trialsize = height(FamilyHistPos.Reconstructed{1,21});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,21}.BadExecution = BadExecution;

%% Participant 22
trialsize = height(FamilyHistPos.Reconstructed{1,22});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,22}.BadExecution = BadExecution;

%% Participant 23
trialsize = height(FamilyHistPos.Reconstructed{1,23});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,23}.BadExecution = BadExecution;

%% Participant 24
trialsize = height(FamilyHistPos.Reconstructed{1,24});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,24}.BadExecution = BadExecution;

%% Participant 25
trialsize = height(FamilyHistPos.Reconstructed{1,25});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(10,1) = 1;

FamilyHistPos.Reconstructed{1,25}.RealReturnAngle(5) = FamilyHistPos.Reconstructed{1,25}.RealReturnAngle(5) + 360;
FamilyHistPos.Reconstructed{1,25}.RealReturnAngle(26) = FamilyHistPos.Reconstructed{1,25}.RealReturnAngle(26) + 360;
FamilyHistPos.Reconstructed{1,25}.RealReturnAngle(36) = FamilyHistPos.Reconstructed{1,25}.RealReturnAngle(36) + 360;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,25}.BadExecution = BadExecution;

%% Participant 26
trialsize = height(FamilyHistPos.Reconstructed{1,26});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,26}.BadExecution = BadExecution;

%% Participant 27
trialsize = height(FamilyHistPos.Reconstructed{1,27});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,27}.RealReturnAngle(23) = FamilyHistPos.Reconstructed{1,27}.RealReturnAngle(23) + 360;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,27}.BadExecution = BadExecution;

%% Participant 28
trialsize = height(FamilyHistPos.Reconstructed{1,28});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,28}.BadExecution = BadExecution;

%% Participant 29
trialsize = height(FamilyHistPos.Reconstructed{1,29});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,29}.RealReturnAngle(15) = FamilyHistPos.Reconstructed{1,29}.RealReturnAngle(15) + 360;
FamilyHistPos.Reconstructed{1,29}.RealReturnAngle(20) = FamilyHistPos.Reconstructed{1,29}.RealReturnAngle(20) + 360;
FamilyHistPos.Reconstructed{1,29}.RealReturnAngle(21) = FamilyHistPos.Reconstructed{1,29}.RealReturnAngle(21) + 360;
FamilyHistPos.Reconstructed{1,29}.RealReturnAngle(23) = FamilyHistPos.Reconstructed{1,29}.RealReturnAngle(23) + 360;
FamilyHistPos.Reconstructed{1,29}.RealReturnAngle(28) = FamilyHistPos.Reconstructed{1,29}.RealReturnAngle(28) + 360;
FamilyHistPos.Reconstructed{1,29}.RealReturnAngle(29) = FamilyHistPos.Reconstructed{1,29}.RealReturnAngle(29) + 360;
FamilyHistPos.Reconstructed{1,29}.RealReturnAngle(32) = FamilyHistPos.Reconstructed{1,29}.RealReturnAngle(32) + 360;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,29}.BadExecution = BadExecution;

%% Participant 30
trialsize = height(FamilyHistPos.Reconstructed{1,30});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(24,1) = 1;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,30}.BadExecution = BadExecution;

%% Participant 31
trialsize = height(FamilyHistPos.Reconstructed{1,31});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,31}.BadExecution = BadExecution;

%% Participant 32
trialsize = height(FamilyHistPos.Reconstructed{1,32});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,32}.RealReturnAngle(7) = FamilyHistPos.Reconstructed{1,32}.RealReturnAngle(7) + 360;
FamilyHistPos.Reconstructed{1,32}.RealReturnAngle(27) = FamilyHistPos.Reconstructed{1,32}.RealReturnAngle(27) + 360;
FamilyHistPos.Reconstructed{1,32}.RealReturnAngle(28) = FamilyHistPos.Reconstructed{1,32}.RealReturnAngle(28) + 360;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,32}.BadExecution = BadExecution;

%% Participant 33
trialsize = height(FamilyHistPos.Reconstructed{1,33});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,33}.BadExecution = BadExecution;

%% Participant 34
trialsize = height(FamilyHistPos.Reconstructed{1,34});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,34}.BadExecution = BadExecution;

%% Participant 35
trialsize = height(FamilyHistPos.Reconstructed{1,35});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,35}.RealReturnAngle(30) = FamilyHistPos.Reconstructed{1,35}.RealReturnAngle(30) + 360;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,35}.BadExecution = BadExecution;

%% Participant 36

trialsize = height(FamilyHistPos.Reconstructed{1,36});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,36}.BadExecution = BadExecution;

%% Participant 37

trialsize = height(FamilyHistPos.Reconstructed{1,37});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,37}.BadExecution = BadExecution;

%% Participant 38

trialsize = height(FamilyHistPos.Reconstructed{1,38});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(2,1) = 1;
BadExecution(5,1) = 1;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,38}.BadExecution = BadExecution;

%% Participant 39

trialsize = height(FamilyHistPos.Reconstructed{1,39});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,39}.BadExecution = BadExecution;

%% Participant 40

trialsize = height(FamilyHistPos.Reconstructed{1,40});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(2,1) = 1;
BadExecution(3,1) = 1;
BadExecution(24,1) = 1;

FamilyHistPos.Reconstructed{1,40}.RealReturnAngle(4) = FamilyHistPos.Reconstructed{1,40}.RealReturnAngle(4) + 360;
FamilyHistPos.Reconstructed{1,40}.RealReturnAngle(21) = FamilyHistPos.Reconstructed{1,40}.RealReturnAngle(21) + 360;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,40}.BadExecution = BadExecution;

%% Participant 41

trialsize = height(FamilyHistPos.Reconstructed{1,41});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,41}.RealReturnAngle(24) = FamilyHistPos.Reconstructed{1,41}.RealReturnAngle(24) + 360;


%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,41}.BadExecution = BadExecution;

%% Participant 42

trialsize = height(FamilyHistPos.Reconstructed{1,42});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,42}.BadExecution = BadExecution;

%% Participant 43

trialsize = height(FamilyHistPos.Reconstructed{1,43});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,43}.BadExecution = BadExecution;

%% Participant 44

trialsize = height(FamilyHistPos.Reconstructed{1,44});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,44}.BadExecution = BadExecution;

%% Participant 45

trialsize = height(FamilyHistPos.Reconstructed{1,45});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,45}.BadExecution = BadExecution;

%% Participant 46

trialsize = height(FamilyHistPos.Reconstructed{1,46});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(21,1) = 1;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,46}.BadExecution = BadExecution;

%% Participant 47

trialsize = height(FamilyHistPos.Reconstructed{1,47});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,47}.BadExecution = BadExecution;

%% Participant 48

trialsize = height(FamilyHistPos.Reconstructed{1,48});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,48}.RealReturnAngle(10) = FamilyHistPos.Reconstructed{1,48}.RealReturnAngle(10) + 360;
FamilyHistPos.Reconstructed{1,48}.RealReturnAngle(28) = FamilyHistPos.Reconstructed{1,48}.RealReturnAngle(28) + 360;


%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,48}.BadExecution = BadExecution;

%% Participant 49

trialsize = height(FamilyHistPos.Reconstructed{1,49});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,49}.BadExecution = BadExecution;

%% Participant 50

trialsize = height(FamilyHistPos.Reconstructed{1,50});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,50}.RealReturnAngle(16) = FamilyHistPos.Reconstructed{1,50}.RealReturnAngle(16) - 360;
FamilyHistPos.Reconstructed{1,50}.RealReturnAngle(36) = FamilyHistPos.Reconstructed{1,50}.RealReturnAngle(36) + 360;


%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,50}.BadExecution = BadExecution;

%% Participant 51

trialsize = height(FamilyHistPos.Reconstructed{1,51});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,51}.BadExecution = BadExecution;

%% Participant 52

trialsize = height(FamilyHistPos.Reconstructed{1,52});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,52}.RealReturnAngle(6) = FamilyHistPos.Reconstructed{1,52}.RealReturnAngle(6) + 360;
FamilyHistPos.Reconstructed{1,52}.RealReturnAngle(12) = FamilyHistPos.Reconstructed{1,52}.RealReturnAngle(12) - 360;
FamilyHistPos.Reconstructed{1,52}.RealReturnAngle(23) = FamilyHistPos.Reconstructed{1,52}.RealReturnAngle(23) + 360;


%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,52}.BadExecution = BadExecution;

%% Participant 53

trialsize = height(FamilyHistPos.Reconstructed{1,53});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(32,1) = 1;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,53}.BadExecution = BadExecution;

%% Participant 54

trialsize = height(FamilyHistPos.Reconstructed{1,54});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,54}.BadExecution = BadExecution;

%% Participant 55

trialsize = height(FamilyHistPos.Reconstructed{1,55});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,55}.BadExecution = BadExecution;

%% Participant 56

trialsize = height(FamilyHistPos.Reconstructed{1,56});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,56}.BadExecution = BadExecution;

%% Participant 57

trialsize = height(FamilyHistPos.Reconstructed{1,57});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,57}.RealReturnAngle(34) = FamilyHistPos.Reconstructed{1,57}.RealReturnAngle(34) + 360;


%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,57}.BadExecution = BadExecution;

%% Participant 58

trialsize = height(FamilyHistPos.Reconstructed{1,58});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,58}.RealReturnAngle(35) = FamilyHistPos.Reconstructed{1,58}.RealReturnAngle(35) + 360;


%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,58}.BadExecution = BadExecution;

%% Participant 59

trialsize = height(FamilyHistPos.Reconstructed{1,59});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,59}.RealReturnAngle(27) = FamilyHistPos.Reconstructed{1,59}.RealReturnAngle(27) + 360;
FamilyHistPos.Reconstructed{1,59}.RealReturnAngle(35) = FamilyHistPos.Reconstructed{1,59}.RealReturnAngle(35) + 360;


%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,59}.BadExecution = BadExecution;

%% Participant 60

trialsize = height(FamilyHistPos.Reconstructed{1,60});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,60}.RealReturnAngle(12) = FamilyHistPos.Reconstructed{1,60}.RealReturnAngle(12) + 360;
FamilyHistPos.Reconstructed{1,60}.RealReturnAngle(32) = FamilyHistPos.Reconstructed{1,60}.RealReturnAngle(32) + 360;


%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,60}.BadExecution = BadExecution;

%% Participant 61

trialsize = height(FamilyHistPos.Reconstructed{1,61});
% Prepare structure
BadExecution = zeros(trialsize,1);

FamilyHistPos.Reconstructed{1,61}.RealReturnAngle(8) = FamilyHistPos.Reconstructed{1,61}.RealReturnAngle(8) + 360;
FamilyHistPos.Reconstructed{1,61}.RealReturnAngle(13) = FamilyHistPos.Reconstructed{1,61}.RealReturnAngle(13) + 360;
FamilyHistPos.Reconstructed{1,61}.RealReturnAngle(27) = FamilyHistPos.Reconstructed{1,61}.RealReturnAngle(27) + 360;
FamilyHistPos.Reconstructed{1,61}.RealReturnAngle(28) = FamilyHistPos.Reconstructed{1,61}.RealReturnAngle(28) + 360;

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,61}.BadExecution = BadExecution;
%% Participant 62

trialsize = height(FamilyHistPos.Reconstructed{1,62});
% Prepare structure
BadExecution = zeros(trialsize,1);

%Save the information in the Reconstructed Structure
FamilyHistPos.Reconstructed{1,62}.BadExecution = BadExecution;

%%
clear trialsize BadExecution