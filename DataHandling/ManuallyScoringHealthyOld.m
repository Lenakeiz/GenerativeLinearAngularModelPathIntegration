%Run Calculate TrackingPath first

%% Participant 1
trialsize = height(HealthyControls.Reconstructed{1,1});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(8,1) = 1;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,1}.BadExecution = BadExecution;

%% Participant 2
trialsize = height(HealthyControls.Reconstructed{1,2});
% Prepare structure
BadExecution = zeros(trialsize,1);

BadExecution(10,1) = 1;

MCINeg.Reconstructed{1,2}.RealReturnAngle(11) = MCINeg.Reconstructed{1,2}.RealReturnAngle(11) + 360;

%Save the information in the Reconstructed Structure
HealthyControls.Reconstructed{1,2}.BadExecution = BadExecution;