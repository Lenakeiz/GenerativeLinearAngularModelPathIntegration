%% Cleaning variables and set intial seed for code reproducibility
clearvars; close all; clc;
rng('default'); 

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
%load('AllDataErrorsPrevent.mat');
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
config.NumParams                                        = 1000;
config.ModelName                                        = "beta_k_g2_g3_sigma_nu";
config.useOoBtrials                                     = true;
config.useTrialFilter                                   = true;

resultfolder = pwd + "/Output/DataFigures/FamilyHistPos_FamilyHistNeg";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% Pos data
FamilyHistPos   = TransformPaths(FamilyHistPos);%transform data
FamilyHistPos   = CalculateTrackingPath(FamilyHistPos, config);
ManuallyScoringFamilyHistPos;
%%
FamilyHistPos.Results = getResultsAllConditions(FamilyHistPos, config);

%% Model fitting for Neg data
FamilyHistNeg   = TransformPaths(FamilyHistNeg);%transform data
FamilyHistNeg   = CalculateTrackingPath(FamilyHistNeg, config);
ManuallyScoringFamilyHistNeg;
%%
FamilyHistNeg.Results = getResultsAllConditions(FamilyHistNeg, config);

%% prepare data for Coco
DataForCoco.FHPos.Code     = FamilyHistPos.Info;
DataForCoco.FHPos.Gender   = FamilyHistPos.Gender;
DataForCoco.FHPos.FlagOoB  = FamilyHistPos.Results.flagOoB;
DataForCoco.FHPos.DistErr  = FamilyHistPos.Results.DistErr;
DataForCoco.FHPos.AngleErr = FamilyHistPos.Results.AngleErr;
DataForCoco.FHPos.ResultsTab   = FamilyHistPos.Reconstructed;
DataForCoco.FHPos.ConditionTab = FamilyHistPos.CondTable;

DataForCoco.FHNeg.Code     = FamilyHistNeg.Info;
DataForCoco.FHNeg.Gender   = FamilyHistNeg.Gender;
DataForCoco.FHNeg.FlagOoB  = FamilyHistNeg.Results.flagOoB;
DataForCoco.FHNeg.DistErr  = FamilyHistNeg.Results.DistErr;
DataForCoco.FHNeg.AngleErr = FamilyHistNeg.Results.AngleErr;
DataForCoco.FHNeg.ResultsTab   = FamilyHistNeg.Reconstructed;
DataForCoco.FHNeg.ConditionTab = FamilyHistNeg.CondTable;

save ./Data/ProcessedDataForCoco.mat DataForCoco

%% Setting colors for using in plots
ColorPattern; 

%% Error plot
config.Gender_Pos = FamilyHistPos.Gender;
config.Gender_Neg = FamilyHistNeg.Gender;
ErrPlot(FamilyHistPos.Results.DistErr, FamilyHistNeg.Results.DistErr, 'dist', config)
ErrPlot(FamilyHistPos.Results.AngleErr, FamilyHistNeg.Results.AngleErr, 'angle', config)
ErrPlot(FamilyHistPos.Results.LocationErr, FamilyHistNeg.Results.LocationErr, 'location', config)

%% Scatter Error Plot
ScatterErrPlot(FamilyHistPos.Results.AngleErr, FamilyHistNeg.Results.AngleErr, FamilyHistPos.Results.flagOoB, FamilyHistNeg.Results.flagOoB, 'male', config)
ScatterErrPlot(FamilyHistPos.Results.AngleErr, FamilyHistNeg.Results.AngleErr, FamilyHistPos.Results.flagOoB, FamilyHistNeg.Results.flagOoB, 'female', config)

%% ThreewayAnova On Data
ThreewayAnova_CocoData(FamilyHistPos.Results.DistErr, FamilyHistNeg.Results.DistErr, config);
ThreewayAnova_CocoData(FamilyHistPos.Results.AngleErr, FamilyHistNeg.Results.AngleErr, config);
ThreewayAnova_CocoData(FamilyHistPos.Results.LocationErr, FamilyHistNeg.Results.LocationErr, config);

%% replicate coc's plots Fig1B PI Error v.s. CAISE DRS
FHPos_DistErrMat = geterrorpercondition(FamilyHistPos.Results.DistErr);
FHNeg_DistErrMat = geterrorpercondition(FamilyHistNeg.Results.DistErr);

FHPos_AngErrMat = geterrorpercondition(FamilyHistPos.Results.AngleErr);
FHNeg_AngErrMat = geterrorpercondition(FamilyHistNeg.Results.AngleErr);

FHPos_LocErrMat = geterrorpercondition(FamilyHistPos.Results.LocationErr);
FHNeg_LocErrMat = geterrorpercondition(FamilyHistNeg.Results.LocationErr);

%get the CAIDE dementia risk score from prevenTab
FHPosDRS = getDRS(FamilyHistPos.Info, preventTab);
FHNegDRS = getDRS(FamilyHistNeg.Info, preventTab);

%%
preventTab.FH; %FH label
preventTab.apoe4; %apoe4 label
nonNanIdx = ~isnan(preventTab.apoe4);

FHTag = preventTab.FH(nonNanIdx);
Apoe4Tag = preventTab.apoe4(nonNanIdx);
DRSTag = preventTab.DRS_noAP(nonNanIdx);  %CAIDE dementia risk score

FHPos_ApoePos = DRSTag(FHTag & Apoe4Tag);   mean(FHPos_ApoePos, 'omitnan')

FHPos_ApoeNeg = DRSTag(FHTag & ~Apoe4Tag);  mean(FHPos_ApoeNeg, 'omitnan')

FHNeg_ApoePos = DRSTag(~FHTag & Apoe4Tag);  mean(FHNeg_ApoePos, 'omitnan')

FHNeg_ApoeNeg = DRSTag(~FHTag & ~Apoe4Tag); mean(FHNeg_ApoeNeg, 'omitnan')

%%
LocErrMat = [FHPos_LocErrMat;FHNeg_LocErrMat];
DRS = [FHPosDRS;FHNegDRS];
y = LocErrMat(:,2)-LocErrMat(:,1);
f = figure('visible','on','Position', [100 100 600 300]);
scatter(DRS, y, 30, "black","filled")
hold on;
h = lsline;
h.LineWidth=2;
h.Color='r';
xlabel("CAIDE Dementia Risk Score")
ylabel("Dist Error (no optical flow-no change)")
%[R,P] = corrcoef(LocErrMat(:,3), DRS, 'Rows','complete');
[R,P] = corrcoef(y, DRS, 'Rows','complete');
title("R="+num2str(R(1,2))+" P="+num2str(P(1,2)))
exportgraphics(f,config.ResultFolder+"/noopticalflow_nochange.png",'Resolution',300);

%% 
AngErrMat = [FHPos_AngErrMat;FHNeg_AngErrMat];
DRS = [FHPosDRS;FHNegDRS];
y = AngErrMat(:,3)-AngErrMat(:,1);
f = figure('visible','on','Position', [100 100 600 300]);
scatter(DRS, y, 30, "black","filled")
hold on;
h = lsline;
h.LineWidth=2;
h.Color='r';
xlabel("CAIDE Dementia Risk Score")
ylabel("Angular Error (no distalcue-no chnange)")
%[R,P] = corrcoef(LocErrMat(:,3), DRS, 'Rows','complete');
[R,P] = corrcoef(y, DRS, 'Rows','complete');
title("R="+num2str(R(1,2))+" P="+num2str(P(1,2)))
exportgraphics(f,config.ResultFolder+"/ang_nodistalcue-nochnage.png",'Resolution',300);

%%
DistErrMat = [FHPos_DistErrMat;FHNeg_DistErrMat];
DRS = [FHPosDRS;FHNegDRS];
y = DistErrMat(:,1);%-DistErrMat(:,1);
f = figure('visible','on','Position', [100 100 600 300]);
scatter(DRS, y, 30, "black","filled")
hold on;
h = lsline;
h.LineWidth=2;
h.Color='r';
xlabel("CAIDE Dementia Risk Score")
ylabel("Dist Error (no chnage)")
[R,P] = corrcoef(LocErrMat(:,1), DRS, 'Rows','complete');
%[R,P] = corrcoef(y, DRS, 'Rows','complete');
title("R="+num2str(R(1,2))+" P="+num2str(P(1,2)))
exportgraphics(f,config.ResultFolder+"/dist_nochnage.png",'Resolution',300);


%% FH+ Apoe+ v.s. FH+ Apoe- v.s. FH- Apoe+ v.s. FH- Apoe- 

% FHPos_ErrMat = FHPos_LocErrMat;
% FHNeg_ErrMat = FHNeg_LocErrMat;

FHPos_ErrMat = FHPos_AngErrMat;
FHNeg_ErrMat = FHNeg_AngErrMat;

% FHPos_ErrMat = FHPos_DistErrMat;
% FHNeg_ErrMat = FHNeg_DistErrMat;


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

nonNanIdx = ~isnan(FHPos_Apoe);
nonNanParams = FHPos_ErrMat(nonNanIdx,:);
nonNanApoe = FHPos_Apoe(nonNanIdx);
nonNanGender = FamilyHistPos.Gender(nonNanIdx);

FHPos_ApoePos_Err = nonNanParams(logical(nonNanApoe),:);
FHPos_ApoeNeg_Err = nonNanParams(~logical(nonNanApoe),:);
FHPos_ApoePos_Gender = nonNanGender(logical(nonNanApoe));
FHPos_ApoeNeg_Gender = nonNanGender(~logical(nonNanApoe));

%extract the subgroups 1,FH+ APoe+ 2,FH+ APoe-....
Gender = FamilyHistNeg.Gender;

nonNanIdx = ~isnan(FHNeg_Apoe);
nonNanParams = FHNeg_ErrMat(nonNanIdx,:);
nonNanApoe = FHNeg_Apoe(nonNanIdx);
nonNanGender = FamilyHistNeg.Gender(nonNanIdx);

FHNeg_ApoePos_Err = nonNanParams(logical(nonNanApoe),Condition);
FHNeg_ApoeNeg_Err = nonNanParams(~logical(nonNanApoe),Condition);
FHNeg_ApoePos_Gender = nonNanGender(logical(nonNanApoe));
FHNeg_ApoeNeg_Gender = nonNanGender(~logical(nonNanApoe));

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
 multicomp_tab123] = ThreewayAnovaOnData_Gender_Condition_FhApoe(FHPos_ApoePos_Err, ...
                                                           FHPos_ApoeNeg_Err, ...
                                                           FHNeg_ApoePos_Err, ...
                                                           FHNeg_ApoeNeg_Err, ...
                                                           config);
%
ErrorBarPlotOfData(FHPos_ApoePos_Err,FHPos_ApoeNeg_Err, FHNeg_ApoePos_Err, FHNeg_ApoeNeg_Err, config)


%% the difference: no distal cue condition - no change condition (baseline)
PosPos_Diff = FHPos_ApoePos_Err(:,3)-FHPos_ApoePos_Err(:,1); %pairwise difference in PosPos 
PosNeg_Diff = FHPos_ApoeNeg_Err(:,3)-FHPos_ApoeNeg_Err(:,1); %pairwise difference in PosNeg 
NegPos_Diff = FHNeg_ApoePos_Err(:,3)-FHNeg_ApoePos_Err(:,1); %pairwise difference in NegPos 
NegNeg_Diff = FHNeg_ApoeNeg_Err(:,3)-FHNeg_ApoeNeg_Err(:,1); %pairwise difference in NegNeg 

TwowayAnovaonDataDiff_Gender_FhApoe(PosPos_Diff,...
                                    PosNeg_Diff,...
                                    NegPos_Diff,...
                                    NegNeg_Diff,...
                                    config);
%
ErrorBarPlotOfDataDiff(PosPos_Diff, PosNeg_Diff, NegPos_Diff, NegNeg_Diff, config);

%%
function ErrorBarPlotOfDataDiff(PosPos_Diff,PosNeg_Diff, NegPos_Diff, NegNeg_Diff, config)

    PosPos_Gender = config.FHPos_ApoePos_Gender;
    PosNeg_Gender = config.FHPos_ApoeNeg_Gender;
    NegPos_Gender = config.FHNeg_ApoePos_Gender;
    NegNeg_Gender = config.FHNeg_ApoeNeg_Gender;
        

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
    PosPos_Male = PosPos_Diff(PosPos_Gender==1); PosPos_Male_mean = mean(PosPos_Male,1,'omitnan'); PosPos_Male_sem = std(PosPos_Male,0,1,'omitnan')/sqrt(size(PosPos_Male,1));
    PosNeg_Male = PosNeg_Diff(PosNeg_Gender==1); PosNeg_Male_mean = mean(PosNeg_Male,1,'omitnan'); PosNeg_Male_sem = std(PosNeg_Male,0,1,'omitnan')/sqrt(size(PosNeg_Male,1));
    NegPos_Male = NegPos_Diff(NegPos_Gender==1); NegPos_Male_mean = mean(NegPos_Male,1,'omitnan'); NegPos_Male_sem = std(NegPos_Male,0,1,'omitnan')/sqrt(size(NegPos_Male,1));
    NegNeg_Male = NegNeg_Diff(NegNeg_Gender==1); NegNeg_Male_mean = mean(NegNeg_Male,1,'omitnan'); NegNeg_Male_sem = std(NegNeg_Male,0,1,'omitnan')/sqrt(size(NegNeg_Male,1));
    
    MaleMean = [NegNeg_Male_mean;NegPos_Male_mean;PosNeg_Male_mean;PosPos_Male_mean];
    MaleSem = [NegNeg_Male_sem;NegPos_Male_sem;PosNeg_Male_sem;PosPos_Male_sem];
    b = bar(1:1:4, MaleMean);

    hold on
    % Plot the errorbars
    errorbar(1:1:4,MaleMean,MaleSem,'k','linestyle','none');
    ylabel("Condition-Baseline")
    title("Male")
    %ylim([-0.2,1.2])
    ylim([-5,20])
    %ylim([-0.2,0.3])

    subplot(1,2,2)
    PosPos_Female = PosPos_Diff(PosPos_Gender==2);PosPos_Female_mean = mean(PosPos_Female,1,'omitnan'); PosPos_Female_sem = std(PosPos_Female,0,1,'omitnan')/sqrt(size(PosPos_Female,1));
    PosNeg_Female = PosNeg_Diff(PosNeg_Gender==2);PosNeg_Female_mean = mean(PosNeg_Female,1,'omitnan'); PosNeg_Female_sem = std(PosNeg_Female,0,1,'omitnan')/sqrt(size(PosNeg_Female,1));
    NegPos_Female = NegPos_Diff(NegPos_Gender==2);NegPos_Female_mean = mean(NegPos_Female,1,'omitnan'); NegPos_Female_sem = std(NegPos_Female,0,1,'omitnan')/sqrt(size(NegPos_Female,1));
    NegNeg_Female = NegNeg_Diff(NegNeg_Gender==2);NegNeg_Female_mean = mean(NegNeg_Female,1,'omitnan'); NegNeg_Female_sem = std(NegNeg_Female,0,1,'omitnan')/sqrt(size(NegNeg_Female,1));

    FemaleMean = [NegNeg_Female_mean;NegPos_Female_mean;PosNeg_Female_mean;PosPos_Female_mean];
    FemaleSem = [NegNeg_Female_sem;NegPos_Female_sem;PosNeg_Female_sem;PosPos_Female_sem];
    b = bar(1:1:4,FemaleMean);

    hold on
    % Plot the errorbars
    errorbar(1:1:4,FemaleMean,FemaleSem,'k','linestyle','none');
    ylabel("Condition-Baseline")
    title("Female")
    %ylim([-0.2,1.2])
    ylim([-5,20])
    %ylim([-0.2,0.3])

    exportgraphics(f,config.ResultFolder+"/DataDiffErrorBar.png",'Resolution',300);
    exportgraphics(f,config.ResultFolder+"/DataDiffErrorBar.pdf",'Resolution',300, 'ContentType','vector');


end

%%
function ErrorBarPlotOfData(PosPos,PosNeg, NegPos, NegNeg, config)

    PosPos_Gender = config.FHPos_ApoePos_Gender;
    PosNeg_Gender = config.FHPos_ApoeNeg_Gender;
    NegPos_Gender = config.FHNeg_ApoePos_Gender;
    NegNeg_Gender = config.FHNeg_ApoeNeg_Gender;

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
    PosPos_Male = PosPos(PosPos_Gender==1,:); PosPos_Male_mean = mean(PosPos_Male,1,'omitnan'); PosPos_Male_sem = std(PosPos_Male,0,1,'omitnan')/sqrt(size(PosPos_Male,1));
    PosNeg_Male = PosNeg(PosNeg_Gender==1,:); PosNeg_Male_mean = mean(PosNeg_Male,1, 'omitnan'); PosNeg_Male_sem = std(PosNeg_Male,0,1,'omitnan')/sqrt(size(PosNeg_Male,1));
    NegPos_Male = NegPos(NegPos_Gender==1,:); NegPos_Male_mean = mean(NegPos_Male,1,'omitnan'); NegPos_Male_sem = std(NegPos_Male,0,1,'omitnan')/sqrt(size(NegPos_Male,1));
    NegNeg_Male = NegNeg(NegNeg_Gender==1,:); NegNeg_Male_mean = mean(NegNeg_Male,1,'omitnan'); NegNeg_Male_sem = std(NegNeg_Male,0,1,'omitnan')/sqrt(size(NegNeg_Male,1));
    
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
    title("Male")
    %ylim([0,2])
    %ylim([0,1.2])

    subplot(1,2,2)
    PosPos_Female = PosPos(PosPos_Gender==2,:);PosPos_Female_mean = mean(PosPos_Female,1,'omitnan'); PosPos_Female_sem = std(PosPos_Female,0,1,'omitnan')/sqrt(size(PosPos_Female,1));
    PosNeg_Female = PosNeg(PosNeg_Gender==2,:);PosNeg_Female_mean = mean(PosNeg_Female,1,'omitnan'); PosNeg_Female_sem = std(PosNeg_Female,0,1,'omitnan')/sqrt(size(PosNeg_Female,1));
    NegPos_Female = NegPos(NegPos_Gender==2,:);NegPos_Female_mean = mean(NegPos_Female,1,'omitnan'); NegPos_Female_sem = std(NegPos_Female,0,1,'omitnan')/sqrt(size(NegPos_Female,1));
    NegNeg_Female = NegNeg(NegNeg_Gender==2,:);NegNeg_Female_mean = mean(NegNeg_Female,1,'omitnan'); NegNeg_Female_sem = std(NegNeg_Female,0,1,'omitnan')/sqrt(size(NegNeg_Female,1));

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
    
    title("Female")
    %ylim([0,2])
    %ylim([0,1.2])
    exportgraphics(f,config.ResultFolder+"/DataErrorBar.png",'Resolution',300);
    exportgraphics(f,config.ResultFolder+"/DataErrorBar.pdf",'Resolution',300, 'ContentType','vector');
end

%%
function Data = addBadExecution(Data)
    %add zero array of BadExecution
    numSubjs = length(Data.Reconstructed);
    for i=1:numSubjs
        trialsize = height(Data.Reconstructed{1,i});
        Data.Reconstructed{1,i}.BadExecution = zeros(trialsize,1);
    end
end

function DRSArray = getDRS(name, preventTab)
    namelist = preventTab.subjectID;
    DRS_noAP = preventTab.DRS_noAP;

    DRSArray = zeros(length(name),1);
    for i=1:length(name)
        name_i = name{i};
        idx = find(strcmp(namelist, name_i));
        DRSArray(i) = DRS_noAP(idx);
    end
end

function ErrMat = geterrorpercondition(Err)
    %swap last two condition
    Condition = [1,3,2];  
    ErrMat = zeros(length(Err{1}),3);
    for i=1:3
        cond = Condition(i);
        Err_cond = Err{cond};
        for id=1:length(Err_cond)
            err = Err_cond{id};
            err = cell2mat(err);
            meanerr = mean(err, 'omitnan');
            ErrMat(id,i) = meanerr;
        end
    end
end

%%
function ScatterErrPlot(FHPosErr, FHNegErr, FHPosFlagOoB, FHNegFlagOoB, gender, config)
    GenderPos = config.Gender_Pos;
    GenderNeg = config.Gender_Neg;

    if gender=="male"
        flag=1;
    else
        flag=2;
    end
    %% Pos 
    PosIdx = GenderPos==flag;

    Pos_InB = zeros(sum(PosIdx),3);
    Pos_OoB = zeros(sum(PosIdx),3);
    Pos_All = zeros(sum(PosIdx),3);

    ScatterPos_OoB = cell(1,3);
    ScatterPos_InB = cell(1,3);
    ScatterPos_All = cell(1,3);

    %%swap last two condition
    Condition = [1,3,2];
    for i=1:3
        cond = Condition(i); 
        ScatterCond_OoB = [];
        ScatterCond_InB = [];
        ScatterCond_All = [];
        Pos_cond_i = FHPosErr{cond}(PosIdx);
        Pos_cond_i_OoB = FHPosFlagOoB{cond}(PosIdx);
        for id=1:length(Pos_cond_i)
            err = Pos_cond_i{id};
            err = cell2mat(err);
            flagOoB = logical(Pos_cond_i_OoB{id});

            InB_err = err(~flagOoB);
            meanerr = mean(InB_err, 'omitnan');
            Pos_InB(id,i) = meanerr;

            OoB_err = err(flagOoB);
            meanerr = mean(OoB_err, 'omitnan');
            Pos_OoB(id,i) = meanerr;

            Pos_All(id,i) = mean(err, 'omitnan');

            ScatterCond_OoB = [ScatterCond_OoB, OoB_err];
            ScatterCond_InB = [ScatterCond_InB, InB_err];
            ScatterCond_All = [ScatterCond_All, err];
        end 

        ScatterPos_OoB{i} = ScatterCond_OoB;
        ScatterPos_InB{i} = ScatterCond_InB;
        ScatterPos_All{i} = ScatterCond_All;
    end

    %% Neg
    NegIdx = GenderNeg==flag;

    Neg_InB = zeros(sum(NegIdx),3);
    Neg_OoB = zeros(sum(PosIdx),3);
    Neg_All = zeros(sum(PosIdx),3);

    ScatterNeg_OoB = cell(1,3);
    ScatterNeg_InB = cell(1,3);
    ScatterNeg_All = cell(1,3);
    
    for i=1:3
        cond = Condition(i);
        ScatterCond_OoB = [];
        ScatterCond_InB = [];
        ScatterCond_All = [];
        Neg_cond_i = FHNegErr{cond}(NegIdx);
        Neg_cond_i_OoB = FHNegFlagOoB{cond}(NegIdx);
        for id=1:length(Neg_cond_i)
            err = Neg_cond_i{id};
            err = cell2mat(err);
            flagOoB = logical(Neg_cond_i_OoB{id});

            InB_err = err(~flagOoB);
            meanerr = mean(InB_err, 'omitnan');
            Neg_InB(id,i) = meanerr;

            OoB_err = err(flagOoB);
            meanerr = mean(OoB_err, 'omitnan');
            Neg_OoB(id,i) = meanerr;

            Neg_All(id,i) = mean(err, 'omitnan');

            ScatterCond_OoB = [ScatterCond_OoB, OoB_err];
            ScatterCond_InB = [ScatterCond_InB, InB_err];
            ScatterCond_All = [ScatterCond_All, err];
        end 

        ScatterNeg_OoB{i} = ScatterCond_OoB;
        ScatterNeg_InB{i} = ScatterCond_InB;
        ScatterNeg_All{i} = ScatterCond_All;
    end

    %% set figure info
    f = figure('visible','on','Position', [100 100 1200 300]);
    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)   

    colorForPos = config.color_scheme_npg(5,:);
    colorForNeg = config.color_scheme_npg(2,:);

    subplot(1,3,1)
    meanPos = mean(Pos_InB, 1, 'omitnan');
    semPos = std(Pos_InB, 1, 'omitnan')/sqrt(size(Pos_InB,1));
    meanNeg = mean(Neg_InB, 1, 'omitnan');
    semNeg = std(Neg_InB, 1, 'omitnan')/sqrt(size(Neg_InB,1));
    errorbar([1,2,3], meanPos, semPos, '-o', 'MarkerSize',10, 'Color', colorForPos, 'LineWidth',3);
    hold on
    errorbar([1,2,3], meanNeg, semNeg, '-o', 'MarkerSize',10, 'Color', colorForNeg, 'LineWidth',3);    

    hold on
    %scatter plot
    for i=1:3
        %Pos InB scatter
        num_points = length(ScatterPos_InB{i});
        x = i*ones(num_points,1)+0.*(rand(num_points,1)-0.5)-0.05; %jitter x
        scatter(x, ScatterPos_InB{i}, 15, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerFaceColor',colorForPos, ...
                'MarkerFaceAlpha',0.2);
        %Neg InB scatter
        num_points = length(ScatterNeg_InB{i});
        x = i*ones(num_points,1)+0.*(rand(num_points,1)-0.5)+0.05; %jitter x
        scatter(x, ScatterNeg_InB{i}, 15, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerFaceColor',colorForNeg, ...
                'MarkerFaceAlpha',0.2); 
    end

    ylabel('Angular error (rads)')
    ylim = [-50,50];
    ytick=[-50,0,50];

    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , (1:3),... 
        'YTick'       , ytick,... 
        'XTickLabel'  , {'No Change','No Optical Flow', 'No Distal Cue'},...
        'XLim'        , [0.5, 3.5],...
        'YLim'        , ylim,...
        'LineWidth'   , .5        );

    title('InB trials')
    legend({'FHPos','FHNeg'})

    subplot(1,3,2)
    meanPos = mean(Pos_OoB, 1, 'omitnan');
    semPos = std(Pos_OoB, 1, 'omitnan')/sqrt(size(Pos_OoB,1));
    meanNeg = mean(Neg_OoB, 1, 'omitnan');
    semNeg = std(Neg_OoB, 1, 'omitnan')/sqrt(size(Neg_OoB,1));
    errorbar([1,2,3], meanPos, semPos, '--o', 'MarkerSize',10, 'Color', colorForPos, 'LineWidth',3);
    hold on
    errorbar([1,2,3], meanNeg, semNeg, '--o', 'MarkerSize',10, 'Color', colorForNeg, 'LineWidth',3);    

    hold on
    %scatter plot
    for i=1:3
        %Pos InB scatter
        num_points = length(ScatterPos_OoB{i});
        x = i*ones(num_points,1)+0.*(rand(num_points,1)-0.5)-0.05; %jitter x
        scatter(x, ScatterPos_OoB{i}, 15, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerFaceColor',colorForPos, ...
                'MarkerFaceAlpha',0.2);
        %Neg InB scatter
        num_points = length(ScatterNeg_OoB{i});
        x = i*ones(num_points,1)+0.*(rand(num_points,1)-0.5)+0.05; %jitter x
        scatter(x, ScatterNeg_OoB{i}, 15, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerFaceColor',colorForNeg, ...
                'MarkerFaceAlpha',0.2); 
    end

    ylabel('Angular error (rads)')
    ylim = [-50,50];
    ytick=[-50,0,50];

    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , (1:3),... 
        'YTick'       , ytick,... 
        'XTickLabel'  , {'No Change','No Optical Flow', 'No Distal Cue'},...
        'XLim'        , [0.5, 3.5],...
        'YLim'        , ylim,...
        'LineWidth'   , .5        );

    title('OoB trials')
    legend({'FHPos','FHNeg'})

    subplot(1,3,3)
    meanPos = mean(Pos_All, 1, 'omitnan');
    semPos = std(Pos_All, 1, 'omitnan')/sqrt(size(Pos_All,1));
    meanNeg = mean(Neg_All, 1, 'omitnan');
    semNeg = std(Neg_All, 1, 'omitnan')/sqrt(size(Neg_All,1));
    errorbar([1,2,3], meanPos, semPos, '--o', 'MarkerSize',10, 'Color', colorForPos, 'LineWidth',3);
    hold on
    errorbar([1,2,3], meanNeg, semNeg, '--o', 'MarkerSize',10, 'Color', colorForNeg, 'LineWidth',3);    

    hold on
    %scatter plot
    for i=1:3
        %Pos InB scatter
        num_points = length(ScatterPos_All{i});
        x = i*ones(num_points,1)+0.*(rand(num_points,1)-0.5)-0.05; %jitter x
        scatter(x, ScatterPos_All{i}, 15, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerFaceColor',colorForPos, ...
                'MarkerFaceAlpha',0.2);
        %Neg InB scatter
        num_points = length(ScatterNeg_All{i});
        x = i*ones(num_points,1)+0.*(rand(num_points,1)-0.5)+0.05; %jitter x
        scatter(x, ScatterNeg_All{i}, 15, ...
                'filled', ...
                'o', ... %marker shape
                'MarkerFaceColor',colorForNeg, ...
                'MarkerFaceAlpha',0.2); 
    end

    ylabel('Angular error (rads)')
    ylim = [-50,50];
    ytick=[-50,0,50];

    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , (1:3),... 
        'YTick'       , ytick,... 
        'XTickLabel'  , {'No Change','No Optical Flow', 'No Distal Cue'},...
        'XLim'        , [0.5, 3.5],...
        'YLim'        , ylim,...
        'LineWidth'   , .5        );

    title('All trials')
    legend({'FHPos','FHNeg'})

end


%%
function ErrPlot(FHPosErr, FHNegErr, type, config)
    %% Pos 
    GenderPos = config.Gender_Pos;
    GenderNeg = config.Gender_Neg;
    PosmaleIdx = GenderPos==1;
    PosfemaleIdx = GenderPos==2;

    PosMale = zeros(sum(PosmaleIdx),3);
    PosFemale = zeros(sum(PosfemaleIdx),3);
    
    %swap last two condition
    Condition = [1,3,2];
    for i=1:3
        cond = Condition(i);
        PosMale_cond_i = FHPosErr{cond}(PosmaleIdx);
        for id=1:length(PosMale_cond_i)
            err = PosMale_cond_i{id};
            err = cell2mat(err);
            meanerr = mean(err, 'omitnan');
            PosMale(id,i) = meanerr;
        end

        PosFemale_cond_i = FHPosErr{cond}(PosfemaleIdx);
        for id=1:length(PosFemale_cond_i)
            err = PosFemale_cond_i{id};
            err = cell2mat(err);
            meanerr = mean(err, 'omitnan');
            PosFemale(id,i) = meanerr;
        end   
    end

%     %swap last two column
%     v = PosMale(:,2); PosMale(:,2)=PosMale(:,3); PosMale(:,3)=v;
%     %swap last two column
%     v = PosFemale(:,2); PosFemale(:,2)=PosFemale(:,3); PosFemale(:,3)=v;
    
    %% Neg
    NegmaleIdx = GenderNeg==1;
    NegfemaleIdx = GenderNeg==2;

    NegMale = zeros(sum(NegmaleIdx),3);
    NegFemale = zeros(sum(NegfemaleIdx),3);

    for i=1:3
        cond = Condition(i);
        NegMale_cond_i = FHNegErr{cond}(NegmaleIdx);
        for id=1:length(NegMale_cond_i)
            err = NegMale_cond_i{id};
            err = cell2mat(err);
            meanerr = mean(err, 'omitnan');
            NegMale(id,i) = meanerr;
        end

        NegFemale_cond_i = FHNegErr{cond}(NegfemaleIdx);
        for id=1:length(NegFemale_cond_i)
            err = NegFemale_cond_i{id};
            err = cell2mat(err);
            meanerr = mean(err, 'omitnan');
            NegFemale(id,i) = meanerr;
        end 
    end

%     %swap last two column
%     v = NegMale(:,2); NegMale(:,2)=NegMale(:,3); NegMale(:,3)=v;
%     %swap last two column
%     v = NegFemale(:,2); NegFemale(:,2)=NegFemale(:,3); NegFemale(:,3)=v;


    %% set figure info
    f = figure('visible','on','Position', [100 100 600 300]);
    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)   

    colorForPos = config.color_scheme_npg(5,:);
    colorForNeg = config.color_scheme_npg(2,:);

    %% male
    subplot(1,2,1);
    meanPosMale = mean(PosMale, 1, 'omitnan');
    semPosMale = std(PosMale, 1, 'omitnan')/sqrt(size(PosMale,1));
    meanNegMale = mean(NegMale, 1, 'omitnan');
    semNegMale = std(NegMale, 1, 'omitnan')/sqrt(size(NegMale,1));
    errorbar([1,2,3], meanPosMale, semPosMale, '-o', 'MarkerSize',10, 'Color', colorForPos, 'LineWidth',3);
    hold on
    errorbar([1,2,3], meanNegMale, semNegMale, '-o', 'MarkerSize',10, 'Color', colorForNeg, 'LineWidth',3);

    if type=="dist"
        ylabel('Distance error (meters)')
        ylim = [0,1.2];
        ytick = [0,0.3,0.6,0.9];
    elseif type=="angle"
        ylabel('Angular error (degrees)')
        %ylim = [-30,30];
        %ytick=[0,10,20];
        ylim = [-20,10];
        ytick=[-20,-10,0,10];
    elseif type=="location"
        ylabel('Location error (meters)')
        ylim = [0,2.0];
        ytick = [0,0.5,1.0,1.5,2.0];        
        %
    end

    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , (1:3),... %'YTick'       , ytick,... 
        'XTickLabel'  , {'No Change','No Optical Flow', 'No Distal Cue'},...
        'XLim'        , [0.5, 3.5],...%'YLim'        , ylim,...
        'LineWidth'   , .5        );

    title('Male')
    legend({'FHPos','FHNeg'})

    %% female
    subplot(1,2,2);
    meanPosFemale = mean(PosFemale, 1, 'omitnan');
    semPosFemale = std(PosFemale, 1, 'omitnan')/sqrt(size(PosFemale,1));
    meanNegFemale = mean(NegFemale, 1, 'omitnan');
    semNegFemale = std(NegFemale, 1, 'omitnan')/sqrt(size(NegFemale,1));
    errorbar([1,2,3], meanPosFemale, semPosFemale, '-o', 'MarkerSize',10, 'Color', colorForPos, 'LineWidth',3);
    hold on
    errorbar([1,2,3], meanNegFemale, semNegFemale, '-o', 'MarkerSize',10, 'Color', colorForNeg, 'LineWidth',3);
  
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , (1:3),... 
        'YTick'       , ytick,... 
        'XTickLabel'  , {'No Change', 'No Optical Flow', 'No Distal Cue'},...
        'XLim'        , [0.5, 3.5],...
        'YLim'        , ylim,...        
        'LineWidth'   , .5        );

    title('Female')
    legend({'FHPos','FHNeg'})

end



%% TwowayAnova and BarScatter Plot
% [anova_tab,multicomp_tab1,~, ~] = TwowayAnova_CocoData(AllFamilyHistPosParams, AllFamilyHistNegParams, config);
% BoxPlotOfFittedParam(AllFamilyHistPosParams, AllFamilyHistNegParams, anova_tab, config);
% BoxPlotOfFittedParamMergeCondition(AllFamilyHistPosParams, AllFamilyHistNegParams, multicomp_tab1, config)

%% ThreewayAnova
% config.Gender_Pos = FamilyHistPos.Gender;
% config.Gender_Neg = FamilyHistNeg.Gender;
% [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab3, multicomp_tab12] = ThreewayAnova_CocoData(AllFamilyHistPosParams, AllFamilyHistNegParams, config);