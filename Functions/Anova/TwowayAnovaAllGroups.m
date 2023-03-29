function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnovaAllGroups(YoungData, HealthyOldData, MCIPosData, MCINegData, MCIUnkData, config)
%% TwowayAnovaOnRealData
% Andrea Castegnaro, UCL, 2022 uceeaca@ucl.ac.uk
% Calculate two way anova between groups of participants
% five groups:         Young / HealthyOld / MCIPos / MCINeg / MCIUnk
% three conditions:    no change / no distal cue /  no optical flow
% Returns anova results and results from multiple comparison tables for
% main effect and interactions
% ===================================================================================

% load output folder from config
resultfolder = config.ResultFolder;
type = config.type;

% Add a subfolder
savefoldername = resultfolder+"/TwowayAnova_"+type+"/";
if ~exist(savefoldername, 'dir')
    mkdir(savefoldername);
end

% processing the data into a long numeric vector
[YoungY, YoungGroupNames, YoungConditionNames]=ReGroupData(YoungData,'Young');
[HealthyOldY, HealthyOldGroupNames, HealthyOldConditionNames]=ReGroupData(HealthyOldData,'HealthyOld');
[MCIPosY, MCIPosGroupNames, MCIPosConditionNames]=ReGroupData(MCIPosData,'MCIPos');
[MCINegY, MCINegGroupNames, MCINegConditionNames]=ReGroupData(MCINegData,'MCIPNeg');
[MCIUnkY, MCIUnkGroupNames, MCIUnkConditionNames]=ReGroupData(MCIUnkData,'MCIPUnk');

    AllY              = [YoungY,             HealthyOldY,             MCIUnkY,             MCINegY,             MCIPosY];
    AllGroupNames     = [YoungGroupNames,    HealthyOldGroupNames,    MCIUnkGroupNames,    MCINegGroupNames,    MCIPosGroupNames];
    AllConditionNames = [YoungConditionNames,HealthyOldConditionNames,MCIUnkConditionNames,MCINegConditionNames,MCIPosConditionNames];
% Do two-way anova with unbalanced design
[p,anova_tab, stats]= anovan(AllY,{AllGroupNames,AllConditionNames},'model','interaction','varnames',{'Groups','Conditions'},'display','on');

% Do multiple comparisons on main effect 1
multicomp_tab1 = multcompare(stats,'Dimension',[1],'CType','bonferroni');
title("Multiple comparisons with bonferroni correction of");
saveas(gcf,savefoldername+"MultiCompME1.png");
close(gcf);

% Do multiple comparisons on main effect 2
multicomp_tab2 = multcompare(stats,'Dimension',[2],'CType','bonferroni');
title("Multiple comparisons with bonferroni correction");
saveas(gcf,savefoldername+"MultiCompME2.png");
close(gcf);

% Do multiple comparisons on main effect 1&2
multicomp_tab12 = multcompare(stats,'Dimension',[1,2],'CType','bonferroni');
title("Multiple comparisons with bonferroni correction");
saveas(gcf,savefoldername+"MultiCompME1ME2.png");
close(gcf);

end