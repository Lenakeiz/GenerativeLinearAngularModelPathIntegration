%% Preparing the data
VAM_PrepareBaseConfig

%% Preprocessing the data
VAM_PreprocessData

%% Setting the model we are interested in
% Eventually modify config paramteters we are interested in. For example
% for this graph we are not interested in running the model with splitted
% conditions so we will set the relative config to false 
% force it to not run
config.useTrialFilter = true;
config.ModelName        =   "beta_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName); %length(config.ParamName); % Set 100 here to avoid producing the model
% Run the model
VAM

%% Model run completed, preparing the data for plotting figures
config.ResultFolder = pwd + "/Output/MRILinearMixedEffectModel/MRIvsModelParams";
% Create storing folder for trajectory if not exist
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

%% Getting Information from results:
YoungControlsParameters   = averageAcrossConditions(YoungControls.Results.estimatedParams);
HealthyControlsParameters = averageAcrossConditions(HealthyControls.Results.estimatedParams);
MCIUnkParameters          = averageAcrossConditions(MCIUnk.Results.estimatedParams);
MCINegParameters          = averageAcrossConditions(MCINeg.Results.estimatedParams);
MCIPosParameters          = averageAcrossConditions(MCIPos.Results.estimatedParams);

%% Preparing linear mixed effect table
HealthyControls_lmetable = createLMETable(HealthyControls.MRI,HealthyControlsParameters);
MCIUnk_lmetable = createLMETable(MCIUnk.MRI,MCIUnkParameters);
MCINeg_lmetable = createLMETable(MCINeg.MRI,MCINegParameters);
MCIPos_lmetable = createLMETable(MCIPos.MRI,MCIPosParameters);

MRIModelParamsDataTable = [HealthyControls_lmetable; MCIUnk_lmetable; MCINeg_lmetable; MCIPos_lmetable];

%Removing nans (either from Nan parameters model or missing mri data)
filter = isnan(MRIModelParamsDataTable.beta) | isnan(MRIModelParamsDataTable.MCI);
MRIModelParamsDataTable = MRIModelParamsDataTable(~filter,:);
clear filter

%% Performing mixedeffectmodel
clc;
rhs = "Sex + Age + norm_ErC + norm_hippocampus + norm_subiculum + norm_isthmuscingulate_volume_Dest + norm_inferiorparietal_volume_Dest + norm_superiorparietal_volume_Dest + (1 | MCI)";

disp("%%%%%%%%%%%%%%% PERFORMING LINEAR MIXED EFFECT MODELS - LARGE MODEL %%%%%%%%%%%%%%%");

lhs = "beta";
[betalme,betastats]   = performLinearMixedEffectModel(MRIModelParamsDataTable, lhs, rhs);
lhs = "g2";
[g2lme,g2stats]       = performLinearMixedEffectModel(MRIModelParamsDataTable, lhs, rhs);
lhs = "g3";
[g3lme,g3stats]       = performLinearMixedEffectModel(MRIModelParamsDataTable, lhs, rhs);
lhs = "sigma";
[sigmalme,sigmastats] = performLinearMixedEffectModel(MRIModelParamsDataTable, lhs, rhs);
lhs = "nu";
[nulme,nustats]       = performLinearMixedEffectModel(MRIModelParamsDataTable, lhs, rhs);

%% Plotting selected variables
close all;

plotInfo.defaultTextSize = 20;
plotInfo.defaultLineSize = 1.3;
plotInfo.titleFontSize = 16;
plotInfo.labelSize = 15;
plotInfo.axisSize = 14;
plotInfo.visible = "off";
plotInfo.ResultFolder = config.ResultFolder;
plotInfo.XLabel = "Subiculum";
plotInfo.YLabel = {'\beta'};
plotSelectedQuantities(MRIModelParamsDataTable.norm_subiculum, MRIModelParamsDataTable.beta, plotInfo);

plotInfo.XLabel = "Hippocampus";
plotInfo.YLabel = '\beta';
plotSelectedQuantities(MRIModelParamsDataTable.norm_hippocampus, MRIModelParamsDataTable.beta, plotInfo);

plotInfo.XLabel = "Inferior parietal";
plotInfo.YLabel = 'g_{2}';
plotSelectedQuantities(MRIModelParamsDataTable.norm_inferiorparietal_volume_Dest, MRIModelParamsDataTable.g2, plotInfo);

plotInfo.XLabel = "Subiculum";
plotInfo.YLabel = 'g_{3}';
plotSelectedQuantities(MRIModelParamsDataTable.norm_subiculum, MRIModelParamsDataTable.g3, plotInfo);

plotInfo.XLabel = "Hippocampus";
plotInfo.YLabel = 'g_{3}';
plotSelectedQuantities(MRIModelParamsDataTable.norm_hippocampus, MRIModelParamsDataTable.g3,plotInfo);

plotInfo.XLabel = "Entorhinal Cortex";
plotInfo.YLabel = '\nu';
plotSelectedQuantities(MRIModelParamsDataTable.norm_ErC, MRIModelParamsDataTable.nu, plotInfo);

plotInfo.XLabel = "Inferior Parietal";
plotInfo.YLabel = '\nu';
plotSelectedQuantities(MRIModelParamsDataTable.norm_inferiorparietal_volume_Dest, MRIModelParamsDataTable.nu, plotInfo);

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
    GroupParamsTable = array2table(GroupParams, "VariableNames", {'beta', 'g2', 'g3', 'sigma', 'nu'});
    dataout = [GroupParamsTable,MRI];
end

%% Perform a linear mized effect model with hypothesis areas
function [lme,stats] = performLinearMixedEffectModel(MRIModelParamsDataTable,leftHandSide, rightHandSide)

    modelFormula = leftHandSide + " ~ " + rightHandSide;
    MRIModelParamsDataTable.CSF = nominal(MRIModelParamsDataTable.CSF);
    MRIModelParamsDataTable.MCI = nominal(MRIModelParamsDataTable.MCI);

    MRIModelParamsDataTable.norm_ErC                          = zscore(MRIModelParamsDataTable.norm_ErC);
    MRIModelParamsDataTable.norm_AntLat                       = zscore(MRIModelParamsDataTable.norm_AntLat);
    MRIModelParamsDataTable.norm_PosMed                       = zscore(MRIModelParamsDataTable.norm_PosMed);
    MRIModelParamsDataTable.norm_entorhinal                   = zscore(MRIModelParamsDataTable.norm_entorhinal);
    MRIModelParamsDataTable.norm_TE35                         = zscore(MRIModelParamsDataTable.norm_TE35);
    MRIModelParamsDataTable.norm_hippocampus                  = zscore(MRIModelParamsDataTable.norm_hippocampus);
    MRIModelParamsDataTable.norm_subiculum                    = zscore(MRIModelParamsDataTable.norm_subiculum);
    MRIModelParamsDataTable.norm_posteriorCingulate           = zscore(MRIModelParamsDataTable.norm_posteriorCingulate);
    MRIModelParamsDataTable.norm_inferiorparietal_volume_Dest = zscore(MRIModelParamsDataTable.norm_inferiorparietal_volume_Dest);
    MRIModelParamsDataTable.norm_isthmuscingulate_volume_Dest = zscore(MRIModelParamsDataTable.norm_isthmuscingulate_volume_Dest);
    MRIModelParamsDataTable.norm_superiorparietal_volume_Dest = zscore(MRIModelParamsDataTable.norm_superiorparietal_volume_Dest);
    MRIModelParamsDataTable.norm_inferiorparietal_volume_DK   = zscore(MRIModelParamsDataTable.norm_inferiorparietal_volume_DK);
    MRIModelParamsDataTable.norm_isthmuscingulate_volume_DK   = zscore(MRIModelParamsDataTable.norm_isthmuscingulate_volume_DK);
    MRIModelParamsDataTable.norm_superiorparietal_volume_DK   = zscore(MRIModelParamsDataTable.norm_superiorparietal_volume_DK);

    disp("%%%%%%%%%%%%%%% LME - " + leftHandSide + " %%%%%%%%%%%%%%%");
    lme = fitglme(MRIModelParamsDataTable,modelFormula)
    stats = anova(lme);

end

%%
function plotSelectedQuantities(x, y, plotInfo)
    % set figure info
    %f = figure('visible','off','Position', [100 100 1000 500]);
    f = figure('visible', plotInfo.visible, 'Position', [100 100 600 600]);
    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)

    mdl = fitlm(x,y);
    plot(mdl);

    xlabel(plotInfo.XLabel, Interpreter="tex");
    ylabel(plotInfo.YLabel, Interpreter="tex");

    title("");

    ax = gca;
    ax.LineWidth = plotInfo.defaultLineSize;
    ax.XLabel.FontSize = plotInfo.labelSize;
    ax.YLabel.FontSize = plotInfo.labelSize;
    ax.XAxis.FontSize = plotInfo.axisSize;
    ax.YAxis.FontSize = plotInfo.axisSize;

    filename = convertCharsToStrings(plotInfo.YLabel) + "vs" + convertCharsToStrings(plotInfo.XLabel);
    
    exportgraphics(f,plotInfo.ResultFolder+"/"+filename+".png",'Resolution',300);
    exportgraphics(f,plotInfo.ResultFolder+"/"+filename+".pdf",'Resolution',300, 'ContentType','vector');

end