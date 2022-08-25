%% Loading NeuroPsych data for all of the older participants filtering for only normalized data
% This script is called from VAM_PrepareBaseConfig.
% 
%
%% Loading MRI data
opts = detectImportOptions("Data\HowettBrain2019_NeuroPsych.csv");
neuroPsychData = readtable("Data\HowettBrain2019_MRI.csv", opts);

clear opts
%% Importing mri data into main structures
Unknown         = LoadMRIIntoGroup(Unknown,mriData);
MCINeg          = LoadMRIIntoGroup(MCINeg,mriData);
MCIPos          = LoadMRIIntoGroup(MCIPos,mriData);
HealthyControls = LoadMRIIntoGroup(HealthyControls,mriData);

clear mriData
%%
function Group = LoadMRIIntoGroup(Group, mriData)
    %% Importing the data
    ids = cellstr(mriData.ID);
    MRI = table(); % not copying the ID
    %MRI.Properties.VariableNames = mriData.Properties.VariableNames(2:end);

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
            MRI = [MRI; mriData(pidx,1:end)];
        else
            emptyTable = array2table(nan(1,width(mriData)));
            emptyTable.Properties.VariableNames = mriData.Properties.VariableNames;
            MRI = [MRI; emptyTable];
        end    
    end

    Group.MRI = MRI;

end
