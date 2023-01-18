%% Preparing the data
VAM_PrepareBaseConfig

%% Preprocessing the data
VAM_PreprocessData

%% Setting the model we are interested in
% Eventually modify config paramteters we are interested in. For example
% for this graph we are not interested in running the model with splitted
% conditions so we will set the relative config to false 
% force it to not run
rng("default")
config.useTrialFilter = true;
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName); %length(config.ParamName); % Set 100 here to avoid producing the model
% Run the model
VAM

%% Model run completed, preparing the data for plotting figures
config.ResultFolder = pwd + "/Output/MRILinearMixedEffectModel/MRIvsModelParams/LinearRegressions/Beta_k_g2_g3_sigma_nu";
% Create storing folder for trajectory if not exist
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

ColorPattern;

%% Getting Information from results:
YoungControlsParameters   = averageAcrossConditions(YoungControls.Results.estimatedParams);
HealthyControlsParameters = averageAcrossConditions(HealthyControls.Results.estimatedParams);
MCIUnkParameters          = averageAcrossConditions(MCIUnk.Results.estimatedParams);
MCINegParameters          = averageAcrossConditions(MCINeg.Results.estimatedParams);
MCIPosParameters          = averageAcrossConditions(MCIPos.Results.estimatedParams);

%% Preparing linear mixed effect table
HealthyControls_lmetable = createLMETable(HealthyControls.MRI,HealthyControlsParameters);
MCIUnk_lmetable          = createLMETable(MCIUnk.MRI,MCIUnkParameters);
MCINeg_lmetable          = createLMETable(MCINeg.MRI,MCINegParameters);
MCIPos_lmetable          = createLMETable(MCIPos.MRI,MCIPosParameters);

MRIModelParamsDataTable  = [HealthyControls_lmetable; MCIUnk_lmetable; MCINeg_lmetable; MCIPos_lmetable];

%Removing nans (either from Nan parameters model or missing mri data)
filter = isnan(MRIModelParamsDataTable.beta) | isnan(MRIModelParamsDataTable.MCI);
MRIModelParamsDataTable = MRIModelParamsDataTable(~filter,:);
clear filter

%% Plotting linear regressions for ROIs
close all;
clc;

plotInfo.defaultTextSize = 14;
plotInfo.defaultLineSize = 1.3;
plotInfo.titleFontSize = 18;
plotInfo.labelSize = 19;
plotInfo.axisSize = 16;
plotInfo.dataSize = 80;
plotInfo.legendFontSize = 15;
plotInfo.visible = "off";
plotInfo.ResultFolder = config.ResultFolder;
plotInfo.color_scheme_group = config.color_scheme_group;

% Getting the indices for the ROI
parameters_label = [{'\beta'}, {'k'}, {'g_2'}, {'g_3'}, {'\sigma'}, {'\nu'}];

mri_indeces = ["norm_ErC"...
    "norm_hippocampus" "norm_superiorparietal_volume_Dest"...
    "norm_inferiorparietal_volume_Dest" "norm_isthmuscingulate_volume_Dest"];
mri_indeces_label = ["Entorhinal Cortex"...
    "Hippocampus" "Superior parietal"...
    "Inferior parietal" "Isthmus cingulate"];

addLine = [0, 1, 1, 1, 0, 0];

for param_i = 1:length(config.ParamName)
    for mri_i = 1:length(mri_indeces)
        plotLinearRegression(...
            MRIModelParamsDataTable,...
            mri_indeces(mri_i),...
            mri_indeces_label(mri_i),...
            config.ParamName(param_i),...
            parameters_label{param_i},...
            addLine(param_i),...
            plotInfo...
        )
    end
end
clear param_i mri_i

clear statuscode sex age

%% get the model parameters, average across the conditions
% remove nans if the row after mergin still contains nans
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

end

%% Creates a unique table between mri and model parameters
function dataout = createLMETable(MRI, GroupParams)
    GroupParamsTable = array2table(GroupParams, "VariableNames", {'beta', 'k', 'g2', 'g3', 'sigma', 'nu'});
    dataout = [GroupParamsTable,MRI];
end

%% Plotting a linear regression using a linear model
function plotLinearRegression(fullTable,x,xLabelPlot,y,yLabelPlot,addLine,plotInfo)

    % set figure info
    f = figure('visible', plotInfo.visible, 'Position', [0 0 500 400]);
    %%% Font type and size setting %%%
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',plotInfo.axisSize)
    set(0,'DefaultTextFontSize',plotInfo.axisSize)

    labels = ["HC", "Unknown", "Negative", "Positive"];
    fullTable.Sex = nominal(fullTable.Sex);
    xData = table2array(fullTable(:,find(strcmp(fullTable.Properties.VariableNames, x))));
    yData = table2array(fullTable(:,find(strcmp(fullTable.Properties.VariableNames, y))));
    csfStatus = fullTable.CSF;

    model_spec = y + " ~ Age + Sex + " + x;
    mdl_full = fitlm(fullTable,model_spec) % To get the pvalue correcting for age and sex
    anova_mdl_full = anova(mdl_full);
    pvalue = anova_mdl_full.pValue(anova_mdl_full.Properties.RowNames == x);
    % Bonferroni corrected
    pvalue = pvalue / 5;
    mdl = fitlm(xData,yData);
    
    max_x = max(xData);
    min_x = min(xData);
    axSc = {};

    hold on;

    for i = 1:4
        axSc{i} = scatter(xData(csfStatus == labels(i)),yData((csfStatus == labels(i))), plotInfo.dataSize);
        axSc{i}.MarkerFaceColor = plotInfo.color_scheme_group(1+i,:);
        axSc{i}.MarkerEdgeColor = plotInfo.color_scheme_group(1+i,:) * 0.8;
        axSc{i}.MarkerFaceAlpha = 0.9; axSc{i}.MarkerEdgeAlpha = 0.9;
    end

    axLine = plot(mdl);
    axLine(1).Marker = 'none';
    axLine(2).Color = [0.2 0.2 0.2];
    axLine(3).Color = [0.2 0.2 0.2];
    axLine(4).Color = [0.2 0.2 0.2];
    axLine(2).LineWidth = plotInfo.defaultLineSize;
    axLine(3).LineWidth = plotInfo.defaultLineSize;
    axLine(4).LineWidth = plotInfo.defaultLineSize;

    xlabel(xLabelPlot, Interpreter="tex");
    ylabel(yLabelPlot, Interpreter="tex");

    title("");
    % pvalue is Bonferroni corrected already
    if(pvalue < 0.001)
        title("***")
    elseif (pvalue < 0.01)
        title("**")
    elseif (pvalue < 0.05)
        title("*")
    end
    legplotInfo = legend([axSc{1}, axSc{2}, axSc{3}, axSc{4}], {'HC' 'MCI unk' 'MCI-' 'MCI+'}, "Location", "northeast", "AutoUpdate", "off");
    legplotInfo.FontSize = plotInfo.legendFontSize;
    
    ax = gca;

    if addLine > 0
        % Plotting y = 1
        tempX = [ax.XLim(1):0.0001:ax.XLim(2)];
        plot(tempX,addLine*ones(length(tempX)),"r--", LineWidth=plotInfo.defaultLineSize-0.2);
    end

    
    ax.LineWidth = plotInfo.defaultLineSize;
    ax.XLabel.FontSize = plotInfo.labelSize;
    ax.YLabel.FontSize = plotInfo.labelSize;
    ax.XAxis.FontSize = plotInfo.axisSize;
    ax.YAxis.FontSize = plotInfo.axisSize;
    ax.XAxis.ExponentMode = "manual";
    ax.XAxis.Exponent = -3;

    hold off;

    filename = convertCharsToStrings(yLabelPlot) + "vs" + convertCharsToStrings(xLabelPlot);
    
    if exist(plotInfo.ResultFolder+"/"+filename+".png","file") == 2
        delete(plotInfo.ResultFolder+"/"+filename+".png");
    end

    exportgraphics(f,plotInfo.ResultFolder+"/"+filename+".png",'Resolution',300);
    exportgraphics(f,plotInfo.ResultFolder+"/"+filename+".pdf",'Resolution',300, 'ContentType','vector');

end