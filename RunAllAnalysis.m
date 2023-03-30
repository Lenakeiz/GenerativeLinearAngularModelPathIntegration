%% Script to run all of the scripts contained in the Main folder. 
% Andrea Castegnaro, UCL, 2022 andrea.castegnaro@ucl.ac.uk
% Be aware this will be time consuming as the code is not optimized.

folder = fullfile(pwd,"Main");
scriptFiles = dir(fullfile(folder, '*.m'));

i_script=1;
while i_script<=numel(scriptFiles)

    scriptName = scriptFiles(i_script).name;
    fprintf('================ Starting execution of %s ================', scriptName);
    
    % Saving locally scriptFiles and i_script counter
    tempFile = fullfile(pwd,"temp.mat");
    save(tempFile,"scriptFiles","i_script");
    
    run(scriptName);
    % Loading scriptFiles and i_script counter
    tempFile = fullfile(pwd,"temp.mat");
    load(tempFile,"scriptFiles","i_script");
    i_script = i_script + 1;

end

tempFile = fullfile(pwd,"temp.mat");
delete(tempFile);