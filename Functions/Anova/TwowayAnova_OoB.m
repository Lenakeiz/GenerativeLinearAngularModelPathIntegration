function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnova_OoB(Young, HO, MCIUnk, MCIPos, MCINeg, config)
%% TwowayAnova_YoungHealthyOldMCICombined
% Zilong Ji, UCL, 2022 zilong.ji@ucl.ac.uk
% Calculate two way anova between groups of participants
% five groups:         Young / HealthyOld / MCIUnk / MCIPos / MCINeg
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

%processing the data into a long numeric vector
[YoungY, YoungGroupNames, YoungConditionNames]=GroupAndRemoveNaN(Young,'Young');
[HOY, HOGroupNames, HOConditionNames]=GroupAndRemoveNaN(HO,'HO');
[MCIUnkY, MCIUnkGroupNames, MCIUnkConditionNames]=GroupAndRemoveNaN(MCIUnk,'MCIUnk');
[MCIPosY, MCIPosGroupNames, MCIPosConditionNames]=GroupAndRemoveNaN(MCIPos,'MCIPos');
[MCINegY, MCINegGroupNames, MCINegConditionNames]=GroupAndRemoveNaN(MCINeg,'MCINeg');


AllY = [YoungY, HOY, MCIUnkY, MCIPosY, MCINegY];
AllGroupNames = [YoungGroupNames,HOGroupNames,MCIUnkGroupNames, MCIPosGroupNames, MCINegGroupNames];
AllConditionNames = [YoungConditionNames,HOConditionNames,MCIUnkConditionNames, MCIPosConditionNames, MCINegConditionNames];

%Do two-way anova with unbalanced design
[p,anova_tab, stats]= anovan(AllY,{AllGroupNames,AllConditionNames},'model','interaction','varnames',{'Groups','Conditions'},'display','on');

%Do multiple comparisons on main effect 1
figure;
multicomp_tab1 = multcompare(stats,'Dimension',[1],'CType','bonferroni');
title("Multiple comparisons with bonferroni correction");
saveas(gcf,savefoldername+"MultiCompME1.png");
%close(gcf);

%Do multiple comparisons on main effect 2
figure;
multicomp_tab2 = multcompare(stats,'Dimension',[2],'CType','bonferroni');
title("Multiple comparisons with bonferroni correction");
saveas(gcf,savefoldername+"MultiCompME2.png");
%close(gcf);

%Do multiple comparisons on main effect 1&2
figure;
multicomp_tab12 = multcompare(stats,'Dimension',[1,2],'CType','bonferroni');
title("Multiple comparisons with bonferroni correction");
saveas(gcf,savefoldername+"MultiCompME1ME2.png");
%close(gcf);

end

%% 
function [Y, GroupNames, ConditionNames]=GroupAndRemoveNaN(Data, groupname)
%% GroupAndRemoveNaN
% Zilong Ji, UCL, 2022, zilong.ji@ucl.ac.uk
% Group all the group data from three environmental conditions into a long
% numeric vector for further two-way anova analysis and also output the
% factor names
% ===================================================================================

Y = [Data(:,1)',Data(:,2)',Data(:,3)'];
GroupNames = [string(repmat({groupname},1,size(Data,1))),...
                 string(repmat({groupname},1,size(Data,1))),...
                 string(repmat({groupname},1,size(Data,1)))];
ConditionNames = [string(repmat({'NoChange'},1,size(Data,1))),...
                 string(repmat({'NoDistalCue'},1,size(Data,1))),...
                 string(repmat({'NoOpticFlow'},1,size(Data,1)))];    


end

