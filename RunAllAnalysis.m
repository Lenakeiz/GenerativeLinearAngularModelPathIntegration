%% Script to run all of the scripts contained in the Main folder. 
% Andrea Castegnaro, UCL, 2022 andrea.castegnaro@ucl.ac.uk
% Be aware this will be time consuming as the code is not optimized.

folder = pwd + "./Main";

scriptFiles = dir(fullfile(folder, '*.m'));

for i=1:numel(scriptFiles)
    scriptName = scriptFiles(i).name;
    run(scriptName);
end

clear folder scriptFiles scriptName