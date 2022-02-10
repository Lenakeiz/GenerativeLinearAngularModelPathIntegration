function TwowayAnovaOnLastTurn(AllYoungTurnAngleDiff, AllHealthyOldTurnAngleDiff, AllMCIPosTurnAngleDiff, AllMCINegTurnAngleDiff, AllMCIUnkTurnAngleDiff, config)
%  Two way ANOVA on turndiff = k3 - mean turn angle

%load configurations necessary for the script
resultfolder = config.ResultFolder;

%create storing folder for trajectory if not exist
savefoldername = resultfolder+"TwowayAnovaOnLastTurn/";
if ~exist(savefoldername, 'dir')
   mkdir(savefoldername);
end
   

%processing the data into a long numeric vector 
[YoungY, YoungGroupNames, YoungConditionNames]=GroupConds(AllYoungTurnAngleDiff,'Young');
[HealthyOldY, HealthyOldGroupNames, HealthyOldConditionNames]=GroupConds(AllHealthyOldTurnAngleDiff,'HealthyOld');
[MCIPosY, MCIPosGroupNames, MCIPosConditionNames]=GroupConds(AllMCIPosTurnAngleDiff,'MCIPos');
[MCINegY, MCINegGroupNames, MCINegConditionNames]=GroupConds(AllMCINegTurnAngleDiff,'MCINeg');
[MCIUnkY, MCIUnkGroupNames, MCIUnkConditionNames]=GroupConds(AllMCIUnkTurnAngleDiff,'MCIUnk');

AllY = [MCIPosY,MCINegY,MCIUnkY,YoungY,HealthyOldY];
AllGroupNames = [MCIPosGroupNames,MCINegGroupNames,MCIUnkGroupNames,YoungGroupNames,HealthyOldGroupNames];
AllConditionNames = [MCIPosConditionNames,MCINegConditionNames,MCIUnkConditionNames,YoungConditionNames,HealthyOldConditionNames];

%Do two-way anova with unbalanced design
[p,tb1, stats]= anovan(AllY,{AllGroupNames,AllConditionNames},'model','interaction','varnames',{'Groups','Conditions'},'display','on');

%Do multiple comparisons on main effect 1
multcompare(stats,'Dimension',[1],'CType','bonferroni');
title("Multiple comparisons with bonferroni correction");
saveas(gcf,savefoldername+"MultiCompME1.png");
close(gcf);

%Do multiple comparisons on main effect 2
multcompare(stats,'Dimension',[2],'CType','bonferroni');
title("Multiple comparisons with bonferroni correction");
saveas(gcf,savefoldername+"MultiCompME2.png");
close(gcf);

%Do multiple comparisons on main effect 1&2
multcompare(stats,'Dimension',[1,2],'CType','bonferroni');
title("Multiple comparisons with bonferroni correction");
saveas(gcf,savefoldername+"MultiCompME1ME2.png");
close(gcf);    

end