%% Cleaning variables and set intial seed for code reproducibility
clearvars; close all; clc;
rng('default'); 

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('AllDataErrorsPreventFH.mat');

%% setting the configuration
config.Speed.alpha                                      = 0.9;    % Paramanter for running speed calculation
config.Speed.timeOffsetAfterFlagReach                   = 1.5;    % Time to track after flag reached in seconds 
config.Speed.smoothWindow                               = 10;     % tracking rate should be 10Hz so 4 secs window is 40 datapoints
config.Speed.velocityCutoff                             = 0.2;    % velocity cutoff to select only the walking part of the reconstructed velocity
config.Speed.timeOffsetForDetectedTemporalWindow        = 0.4;    % time in seconds that will push earlier/ the detected rising edge
config.UseGlobalSearch                                  = true;
config.TrackedInboundAngularDeltaT                      = 1;
config.includeStand                                     = false;
config.useweber                                         = false;  % only true when use weber law in simple generative models
config.Speed.tresholdForBadParticipantL1Recontruction   = 0.0;    % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
config.useOoBtrials                                     = true;
config.useTrialFilter                                   = true;

%% Model fitting
%Model related parameters

% config.ModelName        =   "beta_g2_g3_sigma_nu";
% config.ParamName        =   ["beta", "g2", "g3", "sigma", "nu"];

config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];

config.NumParams        = length(config.ParamName);

resultfolder = pwd + "/Output/ModelFigures/Coco_"+config.ModelName;
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% Model fitting for Pos data
FamilyHistPos   = TransformPaths(FamilyHistPos);%transform data
FamilyHistPos   = CalculateTrackingPath(FamilyHistPos, config);
ManuallyScoringFamilyHistPos;
% fitting 
FamilyHistPos.Results = getResultsAllConditions(FamilyHistPos, config);

%% Model fitting for Neg data
FamilyHistNeg   = TransformPaths(FamilyHistNeg);%transform data
FamilyHistNeg   = CalculateTrackingPath(FamilyHistNeg, config);
ManuallyScoringFamilyHistNeg;
%fitting 
FamilyHistNeg.Results = getResultsAllConditions(FamilyHistNeg, config);

%% extracted the estimated parameters
AllFamilyHistPosParams     =   FamilyHistPos.Results.estimatedParams;
AllFamilyHistNegParams     =   FamilyHistNeg.Results.estimatedParams;

%% Setting colors for using in plots
ColorPattern; 

% %% TwowayAnova and BarScatter Plot
% [anova_tab,multicomp_tab1,~, ~] = TwowayAnova_CocoData(AllFamilyHistPosParams, AllFamilyHistNegParams, config);
% % BoxPlotOfFittedParam(AllFamilyHistPosParams, AllFamilyHistNegParams, anova_tab, config);
% % BoxPlotOfFittedParamMergeCondition(AllFamilyHistPosParams, AllFamilyHistNegParams, multicomp_tab1, config)
% 
% %% ThreewayAnova
% config.Gender_Pos = FamilyHistPos.Gender;
% config.Gender_Neg = FamilyHistNeg.Gender;
% [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab3, multicomp_tab12, multicomp_tab13, multicomp_tab23, multicomp_tab123] = ThreewayAnova_CocoModel(AllFamilyHistPosParams, AllFamilyHistNegParams, config);

%% 4 groups FH+ Apoe+ v.s. FH+ Apoe- v.s. FH- Apoe+ v.s. FH- Apoe-

%find the apoe tag for FH+ and FH-
FHPos_Apoe = zeros(length(FamilyHistPos.Info),1);
for i=1:length(FamilyHistPos.Info)
    name = FamilyHistPos.Info{i};
    idx = find(strcmp(preventTab.subjectID, name));
    %find the apoe tag
    apoe_label = preventTab.apoe4(idx);
    FHPos_Apoe(i) = apoe_label;
end

FHNeg_Apoe = zeros(length(FamilyHistNeg.Info),1);
for i=1:length(FamilyHistNeg.Info)
    name = FamilyHistNeg.Info{i};
    idx = find(strcmp(preventTab.subjectID, name));
    %find the apoe tag
    apoe_label = preventTab.apoe4(idx);
    FHNeg_Apoe(i) = apoe_label;
end

%extract the subgroups 1,FH+ APoe+ 2,FH+ APoe-....
Condition = [1,3,2];  %switch the no distal cue condition with the non optical flow condition
for k=1:3
    cond = Condition(k);
    Params = AllFamilyHistPosParams{1,cond};
    Gender = FamilyHistPos.Gender;

    nonNanIdx = ~isnan(FHPos_Apoe);
    nonNanParams = Params(nonNanIdx,:);
    nonNanApoe = FHPos_Apoe(nonNanIdx);
    nonNanGender = Gender(nonNanIdx);

    FHPos_ApoePos_Param{1,k} = nonNanParams(logical(nonNanApoe),:);
    FHPos_ApoeNeg_Param{1,k} = nonNanParams(~logical(nonNanApoe),:);
    FHPos_ApoePos_Gender = nonNanGender(logical(nonNanApoe));
    FHPos_ApoeNeg_Gender = nonNanGender(~logical(nonNanApoe));
end

%extract the subgroups 1,FH+ APoe+ 2,FH+ APoe-....
for k=1:3
    cond = Condition(k);
    Params = AllFamilyHistNegParams{1,cond};
    Gender = FamilyHistNeg.Gender;

    nonNanIdx = ~isnan(FHNeg_Apoe);
    nonNanParams = Params(nonNanIdx,:);
    nonNanApoe = FHNeg_Apoe(nonNanIdx);
    nonNanGender = Gender(nonNanIdx);

    FHNeg_ApoePos_Param{1,k} = nonNanParams(logical(nonNanApoe),:);
    FHNeg_ApoeNeg_Param{1,k} = nonNanParams(~logical(nonNanApoe),:);
    FHNeg_ApoePos_Gender = nonNanGender(logical(nonNanApoe));
    FHNeg_ApoeNeg_Gender = nonNanGender(~logical(nonNanApoe));
end

config.FHPos_ApoePos_Gender = FHPos_ApoePos_Gender;
config.FHPos_ApoeNeg_Gender = FHPos_ApoeNeg_Gender;
config.FHNeg_ApoePos_Gender = FHNeg_ApoePos_Gender;
config.FHNeg_ApoeNeg_Gender = FHNeg_ApoeNeg_Gender;

% three way anova on gender, codition, and 4 FH+APOE groups
[anova_tab,...
 multicomp_tab1,...
 multicomp_tab2, ...
 multicomp_tab3, ...
 multicomp_tab12, ...
 multicomp_tab13, ...
 multicomp_tab23, ...
 multicomp_tab123] = ThreewayAnova_Gender_Condition_FhApoe(FHPos_ApoePos_Param, ...
                                                           FHPos_ApoeNeg_Param, ...
                                                           FHNeg_ApoePos_Param, ...
                                                           FHNeg_ApoeNeg_Param, ...
                                                           config);
%%
ErrorBarPlotOfFittedParam(FHPos_ApoePos_Param, FHPos_ApoeNeg_Param, FHNeg_ApoePos_Param, FHNeg_ApoeNeg_Param, config);

%% the difference in parmeters no distal cue condition - no change condition (baseline)
PosPos_Param_Diff = FHPos_ApoePos_Param{1,3}-FHPos_ApoePos_Param{1,1}; %pairwise difference in PosPos 
PosNeg_Param_Diff = FHPos_ApoeNeg_Param{1,3}-FHPos_ApoeNeg_Param{1,1}; %pairwise difference in PosNeg 
NegPos_Param_Diff = FHNeg_ApoePos_Param{1,3}-FHNeg_ApoePos_Param{1,1}; %pairwise difference in NegPos 
NegNeg_Param_Diff = FHNeg_ApoeNeg_Param{1,3}-FHNeg_ApoeNeg_Param{1,1}; %pairwise difference in NegNeg 

TwowayAnovaonParameterDiff_Gender_FhApoe(PosPos_Param_Diff,...
                                         PosNeg_Param_Diff,...
                                         NegPos_Param_Diff,...
                                         NegNeg_Param_Diff,...
                                         config);
%
ErrorBarPlotOfParamDiff(PosPos_Param_Diff, PosNeg_Param_Diff, NegPos_Param_Diff, NegNeg_Param_Diff, config);

%%
function ErrorBarPlotOfParamDiff(PosPos_ParamsDiff,PosNeg_ParamsDiff, NegPos_ParamsDiff, NegNeg_ParamsDiff, config)
    ParamName = config.ParamName;

    PosPos_Gender = config.FHPos_ApoePos_Gender;
    PosNeg_Gender = config.FHPos_ApoeNeg_Gender;
    NegPos_Gender = config.FHNeg_ApoePos_Gender;
    NegNeg_Gender = config.FHNeg_ApoeNeg_Gender;

    for ParamIndx=1:length(ParamName)            

        % set figure info
        f = figure('visible','on','Position', [100 100 1000 500]);
        %%% Font type and size setting %%%
        % Using Arial as default because all journals normally require the font to
        % be either Arial or Helvetica
        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)  
        
        subplot(1,2,1)
        PosPos_Male = PosPos_ParamsDiff(PosPos_Gender==1,ParamIndx); PosPos_Male_mean = mean(PosPos_Male,1); PosPos_Male_sem = std(PosPos_Male,0,1)/sqrt(size(PosPos_Male,1));
        PosNeg_Male = PosNeg_ParamsDiff(PosNeg_Gender==1,ParamIndx); PosNeg_Male_mean = mean(PosNeg_Male,1); PosNeg_Male_sem = std(PosNeg_Male,0,1)/sqrt(size(PosNeg_Male,1));
        NegPos_Male = NegPos_ParamsDiff(NegPos_Gender==1,ParamIndx); NegPos_Male_mean = mean(NegPos_Male,1); NegPos_Male_sem = std(NegPos_Male,0,1)/sqrt(size(NegPos_Male,1));
        NegNeg_Male = NegNeg_ParamsDiff(NegNeg_Gender==1,ParamIndx); NegNeg_Male_mean = mean(NegNeg_Male,1); NegNeg_Male_sem = std(NegNeg_Male,0,1)/sqrt(size(NegNeg_Male,1));
        
        MaleMean = [NegNeg_Male_mean;NegPos_Male_mean;PosNeg_Male_mean;PosPos_Male_mean];
        MaleSem = [NegNeg_Male_sem;NegPos_Male_sem;PosNeg_Male_sem;PosPos_Male_sem];
        b = bar(1:1:4, MaleMean);

        hold on
        % Plot the errorbars
        errorbar(1:1:4,MaleMean,MaleSem,'k','linestyle','none');
        title("Male "+ParamName(ParamIndx))

        subplot(1,2,2)
        PosPos_Female = PosPos_ParamsDiff(PosPos_Gender==2,ParamIndx);PosPos_Female_mean = mean(PosPos_Female,1); PosPos_Female_sem = std(PosPos_Female,0,1)/sqrt(size(PosPos_Female,1));
        PosNeg_Female = PosNeg_ParamsDiff(PosNeg_Gender==2,ParamIndx);PosNeg_Female_mean = mean(PosNeg_Female,1); PosNeg_Female_sem = std(PosNeg_Female,0,1)/sqrt(size(PosNeg_Female,1));
        NegPos_Female = NegPos_ParamsDiff(NegPos_Gender==2,ParamIndx);NegPos_Female_mean = mean(NegPos_Female,1); NegPos_Female_sem = std(NegPos_Female,0,1)/sqrt(size(NegPos_Female,1));
        NegNeg_Female = NegNeg_ParamsDiff(NegNeg_Gender==2,ParamIndx);NegNeg_Female_mean = mean(NegNeg_Female,1); NegNeg_Female_sem = std(NegNeg_Female,0,1)/sqrt(size(NegNeg_Female,1));

        FemaleMean = [NegNeg_Female_mean;NegPos_Female_mean;PosNeg_Female_mean;PosPos_Female_mean];
        FemaleSem = [NegNeg_Female_sem;NegPos_Female_sem;PosNeg_Female_sem;PosPos_Female_sem];
        b = bar(1:1:4,FemaleMean);

        hold on
        % Plot the errorbars
        errorbar(1:1:4,FemaleMean,FemaleSem,'k','linestyle','none');
        
        title("Female "+ParamName(ParamIndx))
        exportgraphics(f,config.ResultFolder+"/ParamDiffErrorBar_"+ParamName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/ParamDiffErrorBar_"+ParamName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end


end

%% 
function ErrorBarPlotOfFittedParam(PosPos_Params,PosNeg_Params, NegPos_Params, NegNeg_Params, config)
    numConds = 3; %3 is the condition number
    ParamName = config.ParamName;

    PosPos_Gender = config.FHPos_ApoePos_Gender;
    PosNeg_Gender = config.FHPos_ApoeNeg_Gender;
    NegPos_Gender = config.FHNeg_ApoePos_Gender;
    NegNeg_Gender = config.FHNeg_ApoeNeg_Gender;

    for ParamIndx=1:length(ParamName)
        PosPosParamAllConds = [];
        PosNegParamAllConds = [];
        NegPosParamAllConds = [];
        NegNegParamAllConds = []; %dimension are different, so separate from YoungParams
        
        for TRIAL_FILTER=1:numConds
            %% extract data
            PosPosParam            =   PosPos_Params{TRIAL_FILTER}(:,ParamIndx);
            PosPosParamAllConds    =   [PosPosParamAllConds,PosPosParam];

            PosNegParam            =   PosNeg_Params{TRIAL_FILTER}(:,ParamIndx);
            PosNegParamAllConds    =   [PosNegParamAllConds,PosNegParam];  

            NegPosParam            =   NegPos_Params{TRIAL_FILTER}(:,ParamIndx);
            NegPosParamAllConds    =   [NegPosParamAllConds,NegPosParam];

            NegNegParam            =   NegNeg_Params{TRIAL_FILTER}(:,ParamIndx);
            NegNegParamAllConds    =   [NegNegParamAllConds,NegNegParam];              
        end
        
        %remove Nan rows (nan coz of 1, removing participants with short walking length; 2, not enough trials for parameter estimation)
        PosPosParamAllConds        =   removeNanRows(PosPosParamAllConds);
        PosNegParamAllConds        =   removeNanRows(PosNegParamAllConds);
        NegPosParamAllConds        =   removeNanRows(NegPosParamAllConds);
        NegNegParamAllConds        =   removeNanRows(NegNegParamAllConds);

        % set figure info
        f = figure('visible','on','Position', [100 100 1000 500]);
        %%% Font type and size setting %%%
        % Using Arial as default because all journals normally require the font to
        % be either Arial or Helvetica
        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)  
        
        subplot(1,2,1)
        PosPos_Male = PosPosParamAllConds(PosPos_Gender==1,:); PosPos_Male_mean = mean(PosPos_Male,1); PosPos_Male_sem = std(PosPos_Male,0,1)/sqrt(size(PosPos_Male,1));
        PosNeg_Male = PosNegParamAllConds(PosNeg_Gender==1,:); PosNeg_Male_mean = mean(PosNeg_Male,1); PosNeg_Male_sem = std(PosNeg_Male,0,1)/sqrt(size(PosNeg_Male,1));
        NegPos_Male = NegPosParamAllConds(NegPos_Gender==1,:); NegPos_Male_mean = mean(NegPos_Male,1); NegPos_Male_sem = std(NegPos_Male,0,1)/sqrt(size(NegPos_Male,1));
        NegNeg_Male = NegNegParamAllConds(NegNeg_Gender==1,:); NegNeg_Male_mean = mean(NegNeg_Male,1); NegNeg_Male_sem = std(NegNeg_Male,0,1)/sqrt(size(NegNeg_Male,1));
        
        MaleMean = [NegNeg_Male_mean;NegPos_Male_mean;PosNeg_Male_mean;PosPos_Male_mean];
        MaleSem = [NegNeg_Male_sem;NegPos_Male_sem;PosNeg_Male_sem;PosPos_Male_sem];
        b = bar(MaleMean, 'grouped');

        hold on
        % Calculate the number of groups and number of bars in each group
        [ngroups,nbars] = size(MaleMean);
        % Get the x coordinate of the bars
        x = nan(nbars, ngroups);
        for i = 1:nbars
            x(i,:) = b(i).XEndPoints;
        end
        % Plot the errorbars
        errorbar(x',MaleMean,MaleSem,'k','linestyle','none');
        title(ParamName(ParamIndx)+ " Male")

        subplot(1,2,2)
        PosPos_Female = PosPosParamAllConds(PosPos_Gender==2,:);PosPos_Female_mean = mean(PosPos_Female,1); PosPos_Female_sem = std(PosPos_Female,0,1)/sqrt(size(PosPos_Female,1));
        PosNeg_Female = PosNegParamAllConds(PosNeg_Gender==2,:);PosNeg_Female_mean = mean(PosNeg_Female,1); PosNeg_Female_sem = std(PosNeg_Female,0,1)/sqrt(size(PosNeg_Female,1));
        NegPos_Female = NegPosParamAllConds(NegPos_Gender==2,:);NegPos_Female_mean = mean(NegPos_Female,1); NegPos_Female_sem = std(NegPos_Female,0,1)/sqrt(size(NegPos_Female,1));
        NegNeg_Female = NegNegParamAllConds(NegNeg_Gender==2,:);NegNeg_Female_mean = mean(NegNeg_Female,1); NegNeg_Female_sem = std(NegNeg_Female,0,1)/sqrt(size(NegNeg_Female,1));

        FemaleMean = [NegNeg_Female_mean;NegPos_Female_mean;PosNeg_Female_mean;PosPos_Female_mean];
        FemaleSem = [NegNeg_Female_sem;NegPos_Female_sem;PosNeg_Female_sem;PosPos_Female_sem];
        b = bar(FemaleMean, 'grouped');

        hold on
        % Calculate the number of groups and number of bars in each group
        [ngroups,nbars] = size(FemaleMean);
        % Get the x coordinate of the bars
        x = nan(nbars, ngroups);
        for i = 1:nbars
            x(i,:) = b(i).XEndPoints;
        end
        % Plot the errorbars
        errorbar(x',FemaleMean,FemaleSem,'k','linestyle','none');
        
        title(ParamName(ParamIndx)+ " Female")
        exportgraphics(f,config.ResultFolder+"/ErrorBar_"+ParamName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/ErrorBar_"+ParamName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end


end

%%
function BoxPlotOfFittedParam(AllPosParams, AllNegParams, anova_tab, config)
    
    numConds = 3; %3 is the condition number
    ParamName = config.ParamName;
    for ParamIndx=1:length(ParamName)

        PosParamAllConds = [];
        NegParamAllConds = []; %dimension are different, so separate from YoungParams
        
        for TRIAL_FILTER=1:numConds
            %% extract data
            PosParam            =   AllPosParams{TRIAL_FILTER}(:,ParamIndx);
            PosParamAllConds    =   [PosParamAllConds,PosParam];

            NegParam            =   AllNegParams{TRIAL_FILTER}(:,ParamIndx);
            NegParamAllConds    =   [NegParamAllConds,NegParam];         
        end
        
        %remove Nan rows (nan coz of 1, removing participants with short walking length; 2, not enough trials for parameter estimation)
        PosParamAllConds        =   removeNanRows(PosParamAllConds);
        NegParamAllConds        =   removeNanRows(NegParamAllConds);
    
        %% set figure info
        f = figure('visible','off','Position', [100 100 1000 500]);
        %%% Font type and size setting %%%
        % Using Arial as default because all journals normally require the font to
        % be either Arial or Helvetica
        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)     
        %%% Color definition %%%
        colorForPos = config.color_scheme_npg(3,:);
        colorForNeg = config.color_scheme_npg(5,:);
        
        %set params
        whisker_value = 1.5;
        box_lineWidth = 0.3;
        box_widths_value = 0.2;
        box_color_transparency = 0.5; %faceAlpha
        %center of box (three conditions)
        center_x = [1,2,3];
        shift_value = 0.15; %box shift from center
        median_lineWidth = 2;
        median_color = 'k';
        scatter_jitter_value = 0.1;
        scatter_markerSize=10;
        scatter_marker_edgeColor = 'k';
        scatter_marker_edgeWidth = 0.5;
        scatter_color_transparency = 0.7; %faceAlpha        

        %% boxplot for each column in Pos
        bp1 = boxplot(PosParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x-shift_value);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %% boxplot for each column in Neg
        bp2 = boxplot(NegParamAllConds, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', center_x+shift_value);
        set(bp2,'linewidth',box_lineWidth);

        hold on

        %% Coloring each box
        %findobj first getting the three boxes for Neg (from bp2) from right to left
        %then getting the three boxes for Pos (frm bp1) from right to left
        h = findobj(gca,'Tag','Box'); 
        for i = 1:length(h)
            if i<4  %get the Neg box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForNeg,'FaceAlpha',box_color_transparency);
            else %get the Pos box
                patch(get(h(i),'XData'),get(h(i),'YData'),colorForPos,'FaceAlpha',box_color_transparency);
            end
        end

        %% Adjusting median
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% add scatter plot and the mean of Young
        num_points = size(PosParamAllConds,1);
        for i=1:size(PosParamAllConds,2)
            hold on
            x = i*ones(num_points,1)-shift_value+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, PosParamAllConds(:,i), scatter_markerSize, ...
                    'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForPos, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
            %add errorbar
            mean_Pos = mean(PosParamAllConds(:,i));
            sem_Pos = std(PosParamAllConds(:,i))./sqrt(length(PosParamAllConds(:,i)));
            errorbar(i-shift_value,mean_Pos,sem_Pos,'k','LineStyle','None', 'LineWidth', 2,  'CapSize', 14); 
            hold on
            %add mean point
            scatter(i-shift_value, mean_Pos, 4*scatter_markerSize, 'd',...
                    'filled','MarkerEdgeColor','k', ...
                    'MarkerFaceColor','w', ...
                    'LineWidth',scatter_marker_edgeWidth);
        end

        %% add scatter plot and the mean of Neg
        num_points = size(NegParamAllConds,1);
        for i=1:size(NegParamAllConds,2)
            hold on
            x = i*ones(num_points,1)+shift_value+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
            scatter(x, NegParamAllConds(:,i), scatter_markerSize, ...
                    'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                    'MarkerFaceColor',colorForNeg, ...
                    'MarkerFaceAlpha',scatter_color_transparency,...
                    'LineWidth',scatter_marker_edgeWidth); 
            hold on
            %add errorbar
            mean_Neg = mean(NegParamAllConds(:,i));
            sem_Neg = std(NegParamAllConds(:,i))./sqrt(length(NegParamAllConds(:,i)));
            errorbar(i+shift_value,mean_Neg,sem_Neg,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 14);    
            hold on
            %add mean point
            scatter(i+shift_value, mean_Neg, 4*scatter_markerSize, 'd',...
                    'filled','MarkerEdgeColor','k', ...
                    'MarkerFaceColor','w', ...
                    'LineWidth',scatter_marker_edgeWidth);
        end

        %% Further post-processing the figure

        %calculate te ylim
        alldata = [PosParamAllConds;NegParamAllConds];
        maxdata = max(alldata,[],'all');
        mindata = min(alldata, [], 'all');
        lowupYlim = [mindata-.1*(maxdata-mindata)-eps, maxdata+.1*(maxdata-mindata)+eps]; 
        %eps to make sure when mindata=maxdata, there won't be an error

        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.1 .1 .1], ...
            'YColor'      , [.1 .1 .1], ...
            'XTick'       , (1:3),... 
            'XTickLabel'  , {'No Change','No Distal Cue', 'No Optical Flow'},...
            'XLim'        , [0.5, 3.5],...
            'YLim'        , lowupYlim,... 
            'LineWidth'   , .5        );
        ylabel(ParamName(ParamIndx));
        allpatches = findall(gca,'type','Patch');
        legend(allpatches(1:3:end), {'Pos', 'Neg'}, 'Location','northeast', 'NumColumns',3);

        %extract pvalue for group, conditino and interaction to show on the figure 
        anova_result = anova_tab{ParamIndx};
        group_pvalue = anova_result{2,7};
        condition_pvalue = anova_result{3,7};
        interaction_pvalue = anova_result{4,7};
%         str = {['Group Pvalue = ',sprintf('%.2g',group_pvalue)],...
%                ['Condition Pvalue = ',sprintf('%.2g',condition_pvalue)],...
%                ['Interaction Pvalue = ',sprintf('%.2g',interaction_pvalue)]};
%         annotation('textbox',[0.2 0.6 0.3 0.3],'String',str,'FitBoxToText','on');
        title(strcat(['Group P = ',sprintf('%.2g',group_pvalue)],...
              ['    Condition P = ',sprintf('%.2g',condition_pvalue)],...
              ['    Interaction P = ',sprintf('%.2g',interaction_pvalue)]))

        %% save figure
        exportgraphics(f,config.ResultFolder+"/Box_"+ParamName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/Box_"+ParamName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end

%%
function BoxPlotOfFittedParamMergeCondition(AllPosParams, AllNegParams, multicomp_tab1, config)
    
    numConds = 3; %3 is the condition number
    ParamName = config.ParamName;
    for ParamIndx=1:length(ParamName)

        PosParamAllConds = []; 
        NegParamAllConds = [];

        for TRIAL_FILTER=1:numConds
            %% extract data
            PosParam = AllPosParams{TRIAL_FILTER}(:,ParamIndx);
            PosParamAllConds = [PosParamAllConds,PosParam];

            NegParam = AllNegParams{TRIAL_FILTER}(:,ParamIndx);
            NegParamAllConds = [NegParamAllConds,NegParam];
        end

        %remove Nan rows (nan coz of 1, removing participants with short walking length; 2, not enough trials for parameter estimation)
        PosParamAllConds = removeNanRows(PosParamAllConds);
        NegParamAllConds = removeNanRows(NegParamAllConds);

        PosParamMean = mean(PosParamAllConds, 2);
        NegParamMean = mean(NegParamAllConds, 2);
    
        %% set figure info
        f = figure('visible','off','Position', [100 100 500 500]);
        %%% Font type and size setting %%%
        % Using Arial as default because all journals normally require the font to
        % be either Arial or Helvetica
        set(0,'DefaultAxesFontName','Arial')
        set(0,'DefaultTextFontName','Arial')
        set(0,'DefaultAxesFontSize',12)
        set(0,'DefaultTextFontSize',12)     
        %%% Color definition %%%
        colorForPos = config.color_scheme_npg(6,:);
        colorForNeg = config.color_scheme_npg(3,:);

        %set params
        whisker_value               =   1.5;
        box_lineWidth               =   0.3;
        box_widths_value            =   0.4;
        box_color_transparency      =   0.5; %faceAlpha
        median_lineWidth            =   2;
        median_color                =   'k';
        scatter_jitter_value        =   0.2;
        scatter_markerSize          =   10;
        scatter_marker_edgeColor    =   'k';
        scatter_marker_edgeWidth    =   0.5;
        scatter_color_transparency  =   0.7; %faceAlpha        

        %% boxplot for each column in Pos
        bp1 = boxplot(PosParamMean, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 1);
        set(bp1,'linewidth',box_lineWidth);

        hold on
        %boxplot for each column in Neg
        bp2 = boxplot(NegParamMean, ...
                    'Whisker',whisker_value, ...
                    'symbol','', ... %symbol ='' making outlier invisible
                    'Color','k', ...
                    'Notch','on', ...
                    'widths',box_widths_value,...
                    'positions', 2);
        set(bp2,'linewidth',box_lineWidth);

        %% Coloring each box
        %findobj first getting the box for Neg (from bp2) 
        %then getting the box for Pos(frm bp1)
        h = findobj(gca,'Tag','Box'); 
        %get the Pos box
        patch(get(h(1),'XData'),get(h(1),'YData'),colorForNeg,'FaceAlpha',box_color_transparency);
        %get the Neg box
        patch(get(h(2),'XData'),get(h(2),'YData'),colorForPos,'FaceAlpha',box_color_transparency);

        %% Adjusting median
        h=findobj(gca,'tag','Median');
        for i = 1:length(h)
            h(i).LineWidth = median_lineWidth;
            h(i).Color = median_color;
        end

        %% add scatter plot and the mean of Pos
        num_points = length(PosParamMean);
        hold on
        x = ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, PosParamMean, scatter_markerSize, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForPos, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        %add errorbar
        mean_Pos = mean(PosParamMean);
        sem_Pos= std(PosParamMean)./sqrt(length(PosParamMean));
        errorbar(1,mean_Pos,sem_Pos,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18);    
        hold on
        %add mean point
        scatter(1, mean_Pos, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);

        %% add scatter plot and the mean of Neg
        num_points = size(NegParamMean,1);
        hold on
        x = 2*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, NegParamMean, scatter_markerSize, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',colorForNeg, ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 

        %add errorbar
        mean_Neg = mean(NegParamMean);
        sem_Neg = std(NegParamMean)./sqrt(length(NegParamMean));
        errorbar(2,mean_Neg,sem_Neg,'k','LineStyle','None', 'LineWidth', 2, 'CapSize', 18); 
        hold on
        %add mean point
        scatter(2, mean_Neg, 2*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);        

        %% Further post-processing the figure
        %calculate the ylim
        alldata = [PosParamMean;NegParamMean];
        maxdata = max(alldata,[],'all');
        mindata = min(alldata, [], 'all');
        lowupYlim = [mindata-.1*(maxdata-mindata)-eps, maxdata+.1*(maxdata-mindata)+eps]; 
        %eps to make sure when mindata=maxdata, there won't be an error

        set(gca, ...
            'Box'         , 'off'     , ...
            'TickDir'     , 'out'     , ...
            'TickLength'  , [.01 .01] , ...
            'XColor'      , [.1 .1 .1], ...
            'YColor'      , [.1 .1 .1], ...
            'XTick'       , (1:2),... 
            'XLim'        , [0.5, 2.5],...
            'YLim'        , lowupYlim,...   
            'XTickLabel'  , {'Pos','Neg'},...
            'LineWidth'   , .5        );
        ylabel(ParamName(ParamIndx));

        %extract pvalue for multicomparison of Group effect for showing on the figure
        multicomp_result = multicomp_tab1{ParamIndx};
        Pvalue = multicomp_result(1,6); % Pos vs. Neg see Two-way anova for details
%         str = {['P = ',sprintf('%.2g',Pvalue)]};
%         annotation('textbox',[0.2 0.6 0.3 0.3],'String',str,'FitBoxToText','on');
        title(['P value = ',sprintf('%.2g',Pvalue)])

        %% add sigstar 
        if Pvalue<0.05
            H=sigstar({[1,2]},[Pvalue]);
        end

        %% save figure
        exportgraphics(f,config.ResultFolder+"/ZMergeCondsBox_"+ParamName(ParamIndx)+".png",'Resolution',300);
        exportgraphics(f,config.ResultFolder+"/ZMergeCondsBox_"+ParamName(ParamIndx)+".pdf",'Resolution',300, 'ContentType','vector');

    end
end