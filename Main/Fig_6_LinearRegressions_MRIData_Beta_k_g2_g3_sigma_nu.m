%% Script to create plots from Fig. 6 - Linear regressions models with regions of interest volumes
% Andrea Castegnaro, UCL, 2022 uceeaca@ucl.ac.uk
% Calculate simple linear regressions where the GLAMPI parameters are
% dependent variables and volumes of interest are the independent
% variables. Each linear model includes age and sex as additional
% predictors. We applied Bonferroni correction to the results of the anova
% that tested if the slope of the ROI predictor was different from zero.
% MRI data needs to be loaded for running this script (see details in
% PrepareBaseConfig.

% Preparing the data
GLAMPI_PrepareBaseConfig;

% Preprocessing the data
GLAMPI_PreprocessData;

% Model fitting
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

GLAMPI;

% Preparing output
config.ResultFolder = pwd + "/Output/Fig6";
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Generating color scheme for our paper
ColorPattern; 

% Collecting information from output
YoungControlsParameters   = averageAcrossConditions(YoungControls.Results.estimatedParams);
HealthyControlsParameters = averageAcrossConditions(HealthyControls.Results.estimatedParams);
MCIUnkParameters          = averageAcrossConditions(MCIUnk.Results.estimatedParams);
MCINegParameters          = averageAcrossConditions(MCINeg.Results.estimatedParams);
MCIPosParameters          = averageAcrossConditions(MCIPos.Results.estimatedParams);

% Preparing linear mixed effect table
HealthyControls_lmetable = createLMETable(HealthyControls.MRI,HealthyControlsParameters);
MCIUnk_lmetable          = createLMETable(MCIUnk.MRI,MCIUnkParameters);
MCINeg_lmetable          = createLMETable(MCINeg.MRI,MCINegParameters);
MCIPos_lmetable          = createLMETable(MCIPos.MRI,MCIPosParameters);

MRIModelParamsDataTable  = [HealthyControls_lmetable; MCIUnk_lmetable; MCINeg_lmetable; MCIPos_lmetable];

%Removing nans (either from Nan parameters model or missing mri data)
% Please note this will also account for the exclusion from the
% preprocessing step
filter = isnan(MRIModelParamsDataTable.beta) | isnan(MRIModelParamsDataTable.MCI);
disp("Total Removed samples = " + sum(filter));
MRIModelParamsDataTable = MRIModelParamsDataTable(~filter,:);

disp("Removed helthy controls = " + (height(HealthyControlsParameters) - sum(MRIModelParamsDataTable.CSF == "HC")) ); 
disp("Removed unknown = " + (height(MCIUnkParameters) - sum(MRIModelParamsDataTable.CSF == "Unknown")) ); 
disp("Removed positive = " + (height(MCIPosParameters) - sum(MRIModelParamsDataTable.CSF == "Positive")) ); 
disp("Removed negative = " + (height(MCINegParameters) - sum(MRIModelParamsDataTable.CSF == "Negative")) ); 

clear filter

close all;
clc;

% Parameters set for controlling visual output
plotInfo.defaultTextSize = 12;
plotInfo.defaultLineSize = 1.4;
plotInfo.axisLineSize = 1.2;
plotInfo.LineSizeMdl = 1.5;
plotInfo.titleFontSize = 13;
plotInfo.labelSize = 12;
plotInfo.axisSize = 10;
plotInfo.dataSize = 30;
plotInfo.legendFontSize = 10;
plotInfo.visible = "off";
plotInfo.ResultFolder = config.ResultFolder;
plotInfo.color_scheme_group = config.color_scheme_group;
plotInfo.figurePosition = [200 200 210 210];

parameters_label = [{'\beta'}, {'k'}, {'g_2'}, {'g_3'}, {'\sigma'}, {'\nu'}];

% Labels selected ROI regions, see paper methods for details 
mri_indeces = ["norm_ErC"...
    "norm_hippocampus" "norm_precuneus_volume_Dest"...
    "norm_inferiorparietal_volume_Dest" "norm_isthmuscingulate_volume_Dest"];
mri_indeces_label = ["Entorhinal Cortex"...
    "Hippocampus" "Precuneus"...
    "Inferior parietal" "Isthmus cingulate"];

% Controlling if we should add a reference line when plotting the GLAMPI
% parameter
addLine = [0, 1, 1, 1, 0, 0];

for param_i = 1:length(config.ParamName)    
    for mri_i = 1:length(mri_indeces)

        disp("Fitting " + parameters_label{param_i} + " vs " + mri_indeces_label(mri_i));

        pvalue = fitLinearRegression(...
            MRIModelParamsDataTable,...
            mri_indeces(mri_i),...
            mri_indeces_label(mri_i),...
            config.ParamName(param_i),...
            parameters_label{param_i},...
            addLine(param_i),...
            plotInfo...
        );
        plotLinearRegression(...
            MRIModelParamsDataTable,...
            mri_indeces(mri_i),...
            mri_indeces_label(mri_i),...
            config.ParamName(param_i),...
            parameters_label{param_i},...
            addLine(param_i),...
            plotInfo,...
            pvalue...
        );
    end
end

% Final cleanup to leave workspace as the end of the Preprocessing stage.
% Remove if you want to take a look at the output data.
clearvars -except config YoungControls HealthyControls MCINeg MCIPos MCIUnk

%% ---------------------------------------------------------------------
function dataout = averageAcrossConditions(data)
% get the model parameters, average across the conditions
% remove nans if the row after mergin still contains nans

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

%% ---------------------------------------------------------------------
function dataout = createLMETable(MRI, GroupParams)
% Creates a unique table between mri and model parameters
    GroupParamsTable = array2table(GroupParams, "VariableNames", {'beta', 'k', 'g2', 'g3', 'sigma', 'nu'});
    dataout = [GroupParamsTable,MRI];
end

%%
function pvalue = fitLinearRegression(fullTable,x,xLabelPlot,y,yLabelPlot,addLine,plotInfo)
% Calculated the pvalues from the linear model provided

fullTable.Sex = nominal(fullTable.Sex);
fullTable.CSF = nominal(fullTable.CSF);
fullTable.MCI = nominal(fullTable.MCI);

fullTable.Age                               = zscore(fullTable.Age);  
fullTable.norm_ErC                          = zscore(fullTable.norm_ErC);
fullTable.norm_hippocampus                  = zscore(fullTable.norm_hippocampus);
fullTable.norm_inferiorparietal_volume_Dest = zscore(fullTable.norm_inferiorparietal_volume_Dest);
fullTable.norm_isthmuscingulate_volume_Dest = zscore(fullTable.norm_isthmuscingulate_volume_Dest);
fullTable.norm_precuneus_volume_Dest        = zscore(fullTable.norm_precuneus_volume_Dest);

model_spec = y + " ~ Age + Sex + " + x;
% Let the function display the output
mdl_full = fitlm(fullTable,model_spec)
% Performing anova on coefficients for the slope
anova_mdl_full = anova(mdl_full);

% Output the result for the ROI of interest
pvalue = anova_mdl_full.pValue(anova_mdl_full.Properties.RowNames == x);

% Displaying values that are significant over the Bonferroni correction
if(pvalue < 0.05/5)
    disp("%%%%%%%%%%%%%%%%%%%%%%%");
    disp("%%%%%%%%%%%%%%%%%%%%%%% Survived correction " + xLabelPlot + " " + yLabelPlot + "%%%%%%%%%%%%%%%%%%%%%%%");
    disp("%%%%%%%%%%%%%%%%%%%%%%%");
end

end

%% 
function plotLinearRegression(fullTable,x,xLabelPlot,y,yLabelPlot,addLine,plotInfo,pvalue)
% Plot single linear regression between ROI Volume parameter and GLAMPI
% parameter
    labels = ["HC", "Unknown", "Negative", "Positive"];
    fullTable.Sex = nominal(fullTable.Sex);
    xData = table2array(fullTable(:,find(strcmp(fullTable.Properties.VariableNames, x))));
    yData = table2array(fullTable(:,find(strcmp(fullTable.Properties.VariableNames, y))));
    csfStatus = fullTable.CSF;

    mdl = fitlm(xData,yData);

    % set figure info
    f = figure('visible', plotInfo.visible, 'Position',plotInfo.figurePosition);
    %%% Font type and size setting %%%
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',plotInfo.axisSize)
    set(0,'DefaultTextFontSize',plotInfo.axisSize)

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
    axLine(2).LineWidth = plotInfo.LineSizeMdl;
    axLine(3).LineWidth = plotInfo.LineSizeMdl;
    axLine(4).LineWidth = plotInfo.LineSizeMdl;

    xlabel(xLabelPlot, Interpreter="tex");
    ylabel(yLabelPlot, Interpreter="tex");

    title("");
    
    ax = gca;

    if addLine > 0
        % adding reference line
        tempX = [ax.XLim(1):0.0001:ax.XLim(2)];
        plot(tempX,addLine*ones(length(tempX)),"r--", LineWidth=plotInfo.LineSizeMdl);
    end

    ax.LineWidth = plotInfo.axisLineSize;
    ax.XAxis.FontSize = plotInfo.axisSize;
    ax.YAxis.FontSize = plotInfo.axisSize;
    ax.XLabel.FontSize = plotInfo.labelSize;
    ax.YLabel.FontSize = plotInfo.labelSize;
    ax.XAxis.ExponentMode = "manual";
    ax.XAxis.Exponent = -3;

    hold off;

    legend(ax,"off");

    filename = convertCharsToStrings(yLabelPlot) + "vs" + convertCharsToStrings(xLabelPlot);
    
    exportgraphics(f,plotInfo.ResultFolder+"/"+filename+".png",'Resolution',300);
    exportgraphics(f,plotInfo.ResultFolder+"/"+filename+".pdf",'Resolution',300, 'ContentType','vector');

end