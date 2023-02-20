%% Preparing the data
VAM_PrepareBaseConfig

%% Preprocessing the data
VAM_PreprocessData
 
%% Setting the model we are interested in
% Eventually modify config paramteters we are interested in. For example
% for this graph we are not interested in running the model with splitted
% conditions so we will set the relative config to false 
% force it to not run
rng("default");
config.useTrialFilter = true;
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   100; % Set 100 here to avoid producing the model
% Run the model
VAM

%% Model run completed, preparing the data for plotting figures
config.ResultFolder = pwd + "/Output/Supplementary/SpeedAnalysis";
% Create storing folder for trajectory if not exist
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

%% Genarating color scheme
ColorPattern;

%% Getting speed information
% Averaging across segments L1 and L2 of the outbound path
YoungSpeed = extractSpeedData(YoungControls);
HealthyControlsSpeed = extractSpeedData(HealthyControls);
MCINegSpeed = extractSpeedData(MCINeg);
MCIPosSpeed = extractSpeedData(MCIPos);
MCIUnkSpeed = extractSpeedData(MCIUnk);

% Generate example data
max_element = max([numel(YoungSpeed), numel(HealthyControlsSpeed), numel(MCINegSpeed), numel(MCIPosSpeed), numel(MCIUnkSpeed)]);
YoungSpeed(end+1 : max_element) = nan;
HealthyControlsSpeed(end+1 : max_element) = nan;
MCIUnkSpeed(end+1 : max_element) = nan;
MCINegSpeed(end+1 : max_element) = nan;
MCIPosSpeed(end+1 : max_element) = nan;

SpeedData = [YoungSpeed HealthyControlsSpeed MCIUnkSpeed MCINegSpeed MCIPosSpeed];

figure;

% Create box plot
boxplot(SpeedData,"Notch","on");

% Add data points
hold on;
plot(ones(height(SpeedData),1), SpeedData(:,1), '.', 'MarkerSize', 15, 'Color', 'k');
plot(2*ones(height(SpeedData),1), SpeedData(:,2), '.', 'MarkerSize', 15, 'Color', 'k');
plot(3*ones(height(SpeedData),1), SpeedData(:,3), '.', 'MarkerSize', 15, 'Color', 'k');
plot(4*ones(height(SpeedData),1), SpeedData(:,4), '.', 'MarkerSize', 15, 'Color', 'k');
plot(5*ones(height(SpeedData),1), SpeedData(:,5), '.', 'MarkerSize', 15, 'Color', 'k');
hold off;

ylabel("Speed m/s");
ax = gca;
ax.XAxis.TickLabels = {'Young', 'Healthy Older', 'MCI Unknown', 'MCI Negative', 'MCI Positive'};
title("Outbound path");
%
SpeedDataAnova = [YoungSpeed; HealthyControlsSpeed; MCIUnkSpeed; MCINegSpeed; MCIPosSpeed];

SpeedDataAnovaGroups = [repmat({'Young'}, size(YoungSpeed, 1), 1); ...
          repmat({'Healthy Older'}, size(HealthyControlsSpeed, 1), 1); ...
          repmat({'MCI Unknown'}, size(MCIUnkSpeed, 1), 1); ...
          repmat({'MCI Negative'}, size(MCINegSpeed, 1), 1); ...
          repmat({'MCI Positive'}, size(MCIPosSpeed, 1), 1)];

[p, tbl, stats] = anova1(SpeedDataAnova, SpeedDataAnovaGroups, 'off');
[results, means, ~, ~] = multcompare(stats, 'alpha', 0.05,'Display','off');


%%
function dataout = extractSpeedData(groupData)

    dataout = nan(length(groupData.TrackedL1),2);
    
    for i_participant = 1:length(groupData.TrackedL1)
        % Extracting data calculated from the first segment
        speedTemp = nan(height(groupData.TrackedL1{i_participant}),1);
        for i_trial = 1:height(groupData.TrackedL1{i_participant})
            speedTemp(i_trial) = mean(groupData.TrackedL1{i_participant}{i_trial}.Smoothed_Vel_proj(groupData.TrackedL1{i_participant}{i_trial}.Filtered_Vel_proj),'omitnan');
        end
        dataout(i_participant,1) = mean(speedTemp,'omitnan');

        speedTemp = nan(height(groupData.TrackedL2{i_participant}),1);
        % Extracting data calculated from the second segment
        for i_trial = 1:height(groupData.TrackedL2{i_participant})
            speedTemp(i_trial) = mean(groupData.TrackedL2{i_participant}{i_trial}.Smoothed_Vel_proj(groupData.TrackedL2{i_participant}{i_trial}.Filtered_Vel_proj),'omitnan');
        end
        dataout(i_participant,2) = mean(speedTemp,'omitnan');
    end
    % Calculating average across rows to calculate average speed of participants
    % during the two segments.
    dataout = mean(dataout,2,'omitnan');
end

%% Extracting speed information for each of the walked segment
%
function dataout = averageAcrossConditions(data)

    dataout = [];
    pSize = length(data{1});
    paramsSize = width(data{1});

    for i = 1:pSize
        tempP = [];
        for j = 1:paramsSize
            tempP = [tempP mean([data{1}(i,j) data{2}(i,j) data{3}(i,j)],"omitnan")];
        end
        dataout = [dataout;tempP];
    end

    dataout = removeNanRows(dataout);

end