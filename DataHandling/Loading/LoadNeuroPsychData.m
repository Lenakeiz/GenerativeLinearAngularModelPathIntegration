%% Loading NeuroPsych data for all of the older participants filtering for only normalized data
% This script is called from VAM_PrepareBaseConfig.
% 
%
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

    % Add here the columns you are interested to get from the
    % neuropsych data, if wanted to check for more
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
