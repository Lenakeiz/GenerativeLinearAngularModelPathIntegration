function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnova_Behavioural_MCIOnly(MCIPos, MCINeg, config)
%% TwowayAnova_Behavioural_AllGroups
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
[MCIPosY, MCIPosGroupNames, MCIPosConditionNames]=ReGroupData(MCIPos,'MCIpos');
[MCINegY, MCINegGroupNames, MCINegConditionNames]=ReGroupData(MCINeg,'MCIneg');

AllY              = [MCIPosY,             MCINegY];
AllGroupNames     = [MCIPosGroupNames,    MCINegGroupNames];
AllConditionNames = [MCIPosConditionNames,MCINegConditionNames];
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