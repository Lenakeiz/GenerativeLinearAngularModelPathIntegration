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
config.ResultFolder = pwd + "/Output/MRILinearMixedEffectModel/MRIvsModelParams/Beta_k_g2_g3_sigma_nu";
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
rhs = "Sex + Age + norm_ErC + norm_hippocampus + norm_isthmuscingulate_volume_Dest + norm_inferiorparietal_volume_Dest + norm_precuneus_volume_Dest + (1 | MCI)";

disp("%%%%%%%%%%%%%%% PERFORMING LINEAR MIXED EFFECT MODELS %%%%%%%%%%%%%%%");

lhs = "beta";
[betalme,betastats]   = performLinearMixedEffectModel(MRIModelParamsDataTable, lhs, rhs);
lhs = "k";
[klme,kstats]         = performLinearMixedEffectModel(MRIModelParamsDataTable, lhs, rhs);
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
% plotInfo.XLabel = "Subiculum";
% plotInfo.YLabel = {'\beta'};
% plotSelectedQuantities(MRIModelParamsDataTable.norm_subiculum, MRIModelParamsDataTable.beta, plotInfo);

% plotInfo.XLabel = "Hippocampus";
% plotInfo.YLabel = '\beta';
% plotSelectedQuantities(MRIModelParamsDataTable.norm_hippocampus, MRIModelParamsDataTable.beta, plotInfo);

plotInfo.XLabel = "Inferior parietal";
plotInfo.YLabel = 'g_{2}';
plotSelectedQuantities(MRIModelParamsDataTable.norm_inferiorparietal_volume_Dest, MRIModelParamsDataTable.g2, MRIModelParamsDataTable.CSF, plotInfo, 1);

% plotInfo.XLabel = "Subiculum";
% plotInfo.YLabel = 'g_{3}';
% plotSelectedQuantities(MRIModelParamsDataTable.norm_subiculum, MRIModelParamsDataTable.g3, plotInfo);

plotInfo.XLabel = "Hippocampus";
plotInfo.YLabel = 'g_{3}';
plotSelectedQuantities(MRIModelParamsDataTable.norm_hippocampus, MRIModelParamsDataTable.g3, MRIModelParamsDataTable.CSF, plotInfo, 1);

plotInfo.XLabel = "Entorhinal Cortex";
plotInfo.YLabel = '\nu';
plotSelectedQuantities(MRIModelParamsDataTable.norm_ErC, MRIModelParamsDataTable.nu, MRIModelParamsDataTable.CSF, plotInfo, 0);

plotInfo.XLabel = "Inferior parietal";
plotInfo.YLabel = '\nu';
plotSelectedQuantities(MRIModelParamsDataTable.norm_inferiorparietal_volume_Dest, MRIModelParamsDataTable.nu, MRIModelParamsDataTable.CSF, plotInfo, 0);

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

%% Perform a linear mized effect model with hypothesis areas
function [lme,stats] = performLinearMixedEffectModel(MRIModelParamsDataTable,leftHandSide, rightHandSide)

    modelFormula = leftHandSide + " ~ " + rightHandSide;
    MRIModelParamsDataTable.Sex = nominal(MRIModelParamsDataTable.Sex);
    MRIModelParamsDataTable.CSF = nominal(MRIModelParamsDataTable.CSF);
    MRIModelParamsDataTable.MCI = nominal(MRIModelParamsDataTable.MCI);

    MRIModelParamsDataTable.Age                               = zscore(MRIModelParamsDataTable.Age);  
    MRIModelParamsDataTable.norm_ErC                          = zscore(MRIModelParamsDataTable.norm_ErC);
    MRIModelParamsDataTable.norm_AntLat                       = zscore(MRIModelParamsDataTable.norm_AntLat);
    MRIModelParamsDataTable.norm_PosMed                       = zscore(MRIModelParamsDataTable.norm_PosMed);
    MRIModelParamsDataTable.norm_TE35                         = zscore(MRIModelParamsDataTable.norm_TE35);
    MRIModelParamsDataTable.norm_hippocampus                  = zscore(MRIModelParamsDataTable.norm_hippocampus);
    MRIModelParamsDataTable.norm_subiculum                    = zscore(MRIModelParamsDataTable.norm_subiculum);
    MRIModelParamsDataTable.norm_posteriorCingulate           = zscore(MRIModelParamsDataTable.norm_posteriorCingulate);
    MRIModelParamsDataTable.norm_inferiorparietal_volume_Dest = zscore(MRIModelParamsDataTable.norm_inferiorparietal_volume_Dest);
    MRIModelParamsDataTable.norm_isthmuscingulate_volume_Dest = zscore(MRIModelParamsDataTable.norm_isthmuscingulate_volume_Dest);
    MRIModelParamsDataTable.norm_precuneus_volume_Dest        = zscore(MRIModelParamsDataTable.norm_precuneus_volume_Dest);
    MRIModelParamsDataTable.norm_superiorparietal_volume_DK   = zscore(MRIModelParamsDataTable.norm_superiorparietal_volume_DK);
    MRIModelParamsDataTable.norm_inferiorparietal_volume_DK   = zscore(MRIModelParamsDataTable.norm_inferiorparietal_volume_DK);
    MRIModelParamsDataTable.norm_isthmuscingulate_volume_DK   = zscore(MRIModelParamsDataTable.norm_isthmuscingulate_volume_DK);
    MRIModelParamsDataTable.norm_superiorparietal_volume_DK   = zscore(MRIModelParamsDataTable.norm_superiorparietal_volume_DK);


    disp("%%%%%%%%%%%%%%% LME - " + leftHandSide + " %%%%%%%%%%%%%%%");
    lme = fitlme(MRIModelParamsDataTable,modelFormula)
    disp("%%%%%%%%%%%%%%% Model RSqaured = " + lme.Rsquared.Ordinary + " %%%%%%%%%%%%%%%");    
    stats = anova(lme);

    % Calculating F stat for the model against the null model(random
    % coefficients + intercept + covariates)
    nullModel = leftHandSide + " ~ 1 + (1|MCI)";
    lmeNull = fitlme(MRIModelParamsDataTable,nullModel);
    
    RSS_full = lme.SSE;
    RSS_null = lmeNull.SSE;
    
    n = lme.NumObservations;
    p = lme.NumCoefficients - 2; % only fixed effect (minus intercept)

    dfn = p;
    dfd = n - p - 1;

    F = ((RSS_null - RSS_full) / dfn) / (RSS_full / dfd);

    p_value = 1 - fcdf(F, dfn, dfd);

    [aicNull, bicNull, ~] = aicbic(lmeNull.LogLikelihood,lmeNull.NumPredictors,lmeNull.NumObservations);
    [aicFull, bicFull, ~] = aicbic(lme.LogLikelihood,lme.NumPredictors,lme.NumObservations);

    disp("%%%%%%%%%%%%%%% Fstat(" + dfn + "," + dfd + ") = " + F + " , pValue = " + p_value + " %%%%%%%%%%%%%%%");
    disp("%%%%%%%%%%%%%%% BIC null model = " + bicNull + " ;BIC full model = " + bicFull + " %%%%%%%%%%%%%%%")
end

%%
function plotSelectedQuantities(x, y, csfStatus, plotInfo, addLine)
    % set figure info
    %f = figure('visible','off','Position', [100 100 1000 500]);
    f = figure('visible', plotInfo.visible, 'Position', [0 0 500 400]);
    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',plotInfo.axisSize)
    set(0,'DefaultTextFontSize',plotInfo.axisSize)

    mdl = fitlm(x,y);
    
    labels = ["HC", "Unknown", "Negative", "Positive"];

    hold on;

    for i = 1:4
        axSc{i} = scatter(x(csfStatus == labels(i)),y((csfStatus == labels(i))), plotInfo.dataSize);
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

    xlabel(plotInfo.XLabel, Interpreter="tex");
    ylabel(plotInfo.YLabel, Interpreter="tex");

    title("");
    legplotInfo = legend([axSc{1}, axSc{2}, axSc{3}, axSc{4}], {'HC' 'MCI unk' 'MCI-' 'MCI+'}, "Location", "northeast", "AutoUpdate", "off");
    legplotInfo.FontSize = plotInfo.legendFontSize;
    
    ax = gca;
    ax.LineWidth = plotInfo.defaultLineSize;
    ax.XLabel.FontSize = plotInfo.labelSize;
    ax.YLabel.FontSize = plotInfo.labelSize;
    ax.XAxis.FontSize = plotInfo.axisSize;
    ax.YAxis.FontSize = plotInfo.axisSize;
    ax.XAxis.ExponentMode = "manual";
    ax.XAxis.Exponent = -3;
    ax.XAxis.TickValues = ax.XAxis.TickValues(2:2:end);
    ylim([-0.5 3]);
    ax.YAxis.TickValues = [0 1 2 3];
    ax.YTick = [0 1 2 3];

    if addLine ==  1
        % Plotting y = 1
        tempX = [ax.XAxis.Limits(1):0.0001:ax.XAxis.Limits(2)];
        plot(tempX,ones(length(tempX)),"r--", LineWidth=plotInfo.defaultLineSize-0.2);
    end

    hold off;

    filename = convertCharsToStrings(plotInfo.YLabel) + "vs" + convertCharsToStrings(plotInfo.XLabel);
    
    if exist(plotInfo.ResultFolder+"/"+filename+".png","file") == 2
        delete(plotInfo.ResultFolder+"/"+filename+".png");
    end

    exportgraphics(f,plotInfo.ResultFolder+"/"+filename+".png",'Resolution',300);
    exportgraphics(f,plotInfo.ResultFolder+"/"+filename+".pdf",'Resolution',300, 'ContentType','vector');

end