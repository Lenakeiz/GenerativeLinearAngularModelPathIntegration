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
config.useTrialFilter = false;
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   100; % Set 100 here to avoid producing the model
% Run the model
VAM

%% Model run completed, preparing the data for plotting figures
config.ResultFolder = pwd + "/Output/Supplementary/TimeAnalysis";
% Create storing folder for trajectory if not exist
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

%% Genarating color scheme
ColorPattern;

%% Getting Information from results:
YoungControlsParameters   = averageAcrossConditions(YoungControls.Results.estimatedParams);
HealthyControlsParameters = averageAcrossConditions(HealthyControls.Results.estimatedParams);
MCIUnkParameters          = averageAcrossConditions(MCIUnk.Results.estimatedParams);
MCINegParameters          = averageAcrossConditions(MCINeg.Results.estimatedParams);
MCIPosParameters          = averageAcrossConditions(MCIPos.Results.estimatedParams);

%% Getting speed information
% Averaging across segments L1 and L2 of the outbound path
YoungOutboundDuration = extractTimeData(YoungControls);
HealthyControlsOutboundDuration = extractTimeData(HealthyControls);
MCINegOutboundDuration = extractTimeData(MCINeg);
MCIPosOutboundDuration = extractTimeData(MCIPos);
MCIUnkOutboundDuration = extractTimeData(MCIUnk);
%%
% Generate example data
max_element = max([numel(YoungOutboundDuration),...
    numel(HealthyControlsOutboundDuration),...
    numel(MCINegOutboundDuration),...
    numel(MCIPosOutboundDuration),...
    numel(MCIUnkOutboundDuration)]);
YoungOutboundDuration(end+1 : max_element) = nan;
HealthyControlsOutboundDuration(end+1 : max_element) = nan;
MCIUnkOutboundDuration(end+1 : max_element) = nan;
MCINegOutboundDuration(end+1 : max_element) = nan;
MCIPosOutboundDuration(end+1 : max_element) = nan;

OutboundDurationData = [YoungOutboundDuration HealthyControlsOutboundDuration MCIUnkOutboundDuration MCINegOutboundDuration MCIPosOutboundDuration];

% Create box plot
boxplot(OutboundDurationData);

% Add data points
hold on;
plot(ones(height(OutboundDurationData),1), OutboundDurationData(:,1), '.', 'MarkerSize', 15, 'Color', 'k');
plot(2*ones(height(OutboundDurationData),1), OutboundDurationData(:,2), '.', 'MarkerSize', 15, 'Color', 'k');
plot(3*ones(height(OutboundDurationData),1), OutboundDurationData(:,3), '.', 'MarkerSize', 15, 'Color', 'k');
plot(4*ones(height(OutboundDurationData),1), OutboundDurationData(:,4), '.', 'MarkerSize', 15, 'Color', 'k');
plot(5*ones(height(OutboundDurationData),1), OutboundDurationData(:,5), '.', 'MarkerSize', 15, 'Color', 'k');
hold off;

ylabel("Duration (s)");
ax = gca;
ax.XAxis.TickLabels = {'Young', 'Healthy Older', 'MCI Unknown', 'MCI Negative', 'MCI Positive'};
title("Outbound path");
%
OutboundDurationDataAnova = [YoungOutboundDuration; HealthyControlsOutboundDuration; MCIUnkOutboundDuration; MCINegOutboundDuration; MCIPosOutboundDuration];

OutboundDurationDataAnovaGroups = [repmat({'Young'}, size(YoungOutboundDuration, 1), 1); ...
          repmat({'Healthy Older'}, size(HealthyControlsOutboundDuration, 1), 1); ...
          repmat({'MCI Unknown'}, size(MCIUnkOutboundDuration, 1), 1); ...
          repmat({'MCI Negative'}, size(MCINegOutboundDuration, 1), 1); ...
          repmat({'MCI Positive'}, size(MCIPosOutboundDuration, 1), 1)];

[p, tbl, stats] = anova1(OutboundDurationDataAnova, OutboundDurationDataAnovaGroups, 'off');
[results, means, ~, ~] = multcompare(stats, 'alpha', 0.05,'Display','off');


%%
function dataout = extractTimeData(groupData)

    dataout = [];
    
    for i_participant = 1:length(groupData.Results.L1Dur)
        if(~isempty(groupData.Results.L1Dur{i_participant}))
            l1_duration = cell2mat(groupData.Results.L1Dur{i_participant});
            l2_duration = cell2mat(groupData.Results.L2Dur{i_participant});
            dataout = [dataout; mean(l1_duration+l2_duration,'omitnan')];
        end
    end
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