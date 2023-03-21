%% DescriptiveStatitisticsPreprocessing
% This script will output the number of used trials at the end of the
% preprocessing stage. It will also output the number of outbound trials
% proportion for each separate group

%% Preparing the data
VAM_PrepareBaseConfig

%% Preprocessing the data
VAM_PreprocessData
 
%% Preparing and fitting the model.
% For more information see  
rng("default");
config.useTrialFilter = true;
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName); % Set 100 here to avoid producing the model
% Run the model
VAM

%% Clearing the workspace
OutputPreprocessingDescriptive(YoungControls,"Young controls");
OutputPreprocessingDescriptive(HealthyControls,"Healthy Controls");
OutputPreprocessingDescriptive(MCIUnk,"MCI Unk");
OutputPreprocessingDescriptive(MCINeg,"MCI Neg");
OutputPreprocessingDescriptive(MCIPos,"MCI pos");
clearvars -except config YoungControls HealthyControls MCIUnk MCINeg MCIPos

%% Takes a sigle group and calcutate the descriptive statistics out of it
function OutputPreprocessingDescriptive(group, groupName)
% Looping through each of the groups to calculate properties
n_participant = length(group.Reconstructed);
excludedTrials = 0;
outOfBoundTrials = 0;
totalTrials = 0;
% We can take the following information from any participant
experimentTrials = length(group.FlagPos{1});

for i=1:n_participant
    totalTrials = totalTrials + experimentTrials;

    if(ismember(i,group.BadPptIdxs))
        excludedTrials = excludedTrials + experimentTrials;
        continue;
    end
    excludedTrials = excludedTrials + sum(group.Reconstructed{i}.BadExecution);
    % We count out of bound trials only the one that are not also excluded.
    outOfBoundTrials = outOfBoundTrials + sum(group.CondTable{i}.OutOfBound - group.Reconstructed{i}.BadExecution > 0);
end

ratioExcluded = (excludedTrials/totalTrials) * 100;
ratioOob = (outOfBoundTrials/totalTrials) * 100;

disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
disp(sprintf('%%%%%%%%%%%%%% Group Name: %s %%%%%%%%%%%%%%%%\nTotal Trials = %d\nExcluded Trials = %d (%.2f%%)\nOut of bound Trials = %d (%.2f%%)', ...
groupName, totalTrials, excludedTrials, ratioExcluded, outOfBoundTrials, ratioOob))
end




