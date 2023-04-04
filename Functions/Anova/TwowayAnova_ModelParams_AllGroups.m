function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnova_ModelParams_AllGroups(YoungData, HealthyOldData, MCIUnkData, MCINegData, MCIPosData, config)
%% TwowayAnova_ModelParams_AllGroups
% Andrea Castegnaro, UCL, 2023 uceeaca@ucl.ac.uk
% Calculate two way anova between groups of participants
% five groups:         Young / HealthyOld / MCIunk / MCINeg MCIPos
% three conditions:    no change / no distal cue /  no optical flow
% Expects a set of parameters for each group as an input
% Returns anova results and results from multiple comparison tables for
% main effect and interactions.
% ===================================================================================

%load output
resultfolder = config.ResultFolder;

%add subfolder
savefoldername = resultfolder+"/TwowayAnova/";
if ~exist(savefoldername, 'dir')
    mkdir(savefoldername);
end

ParamName = config.ParamName;
param_nums = length(ParamName);

anova_tab = cell(0);
multicomp_tab1 = cell(0);
multicomp_tab2 = cell(0);
multicomp_tab12 = cell(0);

%for each parameter, calculate the anova and multiple test
for param_idx=1:param_nums

    param_name = ParamName(param_idx);

    %processing the data into a long numeric vector
    [YoungY, YoungGroupNames, YoungConditionNames]=GroupAndRemoveNaN(YoungData,param_idx,'Young');
    [HealthyOldY, HealthyOldGroupNames, HealthyOldConditionNames]=GroupAndRemoveNaN(HealthyOldData,param_idx,'HealthyOld');
    [MCIUnkY, MCIUnkGroupNames, MCIUnkConditionNames]=GroupAndRemoveNaN(MCIUnkData,param_idx,'Unknown');
    [MCINegY, MCINegGroupNames, MCINegConditionNames]=GroupAndRemoveNaN(MCINegData,param_idx,'Negative');
    [MCIPosY, MCIPosGroupNames, MCIPosConditionNames]=GroupAndRemoveNaN(MCIPosData,param_idx,'Positive');

    AllY = [YoungY, HealthyOldY, MCIUnkY, MCINegY, MCIPosY];
    AllGroupNames = [YoungGroupNames, HealthyOldGroupNames, MCIUnkGroupNames, MCINegGroupNames, MCIPosGroupNames];
    AllConditionNames = [YoungConditionNames, HealthyOldConditionNames, MCIUnkConditionNames, MCINegConditionNames, MCIPosConditionNames];

    %Do two-way anova with unbalanced design
    [p,tb1, stats]= anovan(AllY,{AllGroupNames,AllConditionNames},'model','interaction','varnames',{'Groups','Conditions'},'display','on');
    anova_tab{param_idx} = tb1;

    %Do multiple comparisons on main effect 1
    result = multcompare(stats,'Dimension',[1],'CType','bonferroni');
    multicomp_tab1{param_idx} = result;
    title("Multiple comparisons with bonferroni correction of : "+param_name);
    saveas(gcf,savefoldername+"MultiCompME1_"+param_name+".png");
    close(gcf);

    %Do multiple comparisons on main effect 2
    result = multcompare(stats,'Dimension',[2],'CType','bonferroni');
    multicomp_tab2{param_idx} = result;
    title("Multiple comparisons with bonferroni correction of parameter: "+param_name);
    saveas(gcf,savefoldername+"MultiCompME2_"+param_name+".png");
    close(gcf);

    %Do multiple comparisons on main effect 1&2
    result = multcompare(stats,'Dimension',[1,2],'CType','bonferroni');
    multicomp_tab12{param_idx} = result;
    title("Multiple comparisons with bonferroni correction of parameter: "+param_name);
    saveas(gcf,savefoldername+"MultiCompME1ME2_"+param_name+".png");
    close(gcf);
end

end

