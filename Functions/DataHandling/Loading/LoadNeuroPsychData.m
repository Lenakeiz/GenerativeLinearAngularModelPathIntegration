%% Loading neuropsycholigcal data for elderly participants
% Andrea Castegnaro, UCL, 2022, uceeaca@ucl.ac.uk
% Loads neuropsychological tests acquitred in Howett et al., 2019 study.
% Tries to link the generated identifier from the extracted data with the
% one from the behavioural data.
% The data extracted here is not currently used in the present study.
% This script is called from VAM_PrepareBaseConfig.

%% Loading MRI data
opts = detectImportOptions("./Data/HowettBrain2019_NeuroPsych.csv");
neuroPsychData = readtable("./Data/HowettBrain2019_NeuroPsych.csv", opts);

clear opts
%% Importing mri data into main structures
Unknown         = AddNeuroPsychData(Unknown,neuroPsychData);
MCINeg          = AddNeuroPsychData(MCINeg,neuroPsychData);
MCIPos          = AddNeuroPsychData(MCIPos,neuroPsychData);
HealthyControls = AddNeuroPsychData(HealthyControls,neuroPsychData);

clear neuroPsychData
%%
function Group = AddNeuroPsychData(Group, neuroPsychData)
    %% Importing the data
    ids = cellstr(neuroPsychData.ID);
    NeuroPsych = table();

    % Extracting specific tests indexed by column within the dataset
    neuroPsychIndices = [15 16 20];

    for i = 1:length(Group.Info)
        extractID = Group.Info{i};
        extractID = extractAfter(extractID,"A093863_");
        pid       = extractBefore(extractID,"_");
        pidInt = str2num(pid);

        searchpid = "";
        if(pidInt < 10)
            searchpid = append("A093863_0",pid);
        else
            searchpid = append("A093863_",pid);
        end

        pidx = find(strcmp(ids,searchpid));
        if(~isempty(pidx))

            % collected neuropsych tests
            NeuroPsych = [NeuroPsych; neuroPsychData(pidx,neuroPsychIndices)];
        else
            emptyTable = array2table(nan(1,length(neuroPsychIndices)));
            emptyTable.Properties.VariableNames = neuroPsychData.Properties.VariableNames(neuroPsychIndices);
            NeuroPsych = [NeuroPsych; emptyTable];
        end    
    end

    Group.NeuroPsych = NeuroPsych;

end
