%% Two-way anova on real data: distance, angle or mean walking speed
% five groups:         Young / HealthyOld / MCIPos / MCINeg / MCIUnk
% three conditions:    no change / no distal cue /  no optical flow
function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnovaOnRealData(YoungData, HealthyOldData, MCIPosData, MCINegData, MCIUnkData, config)

    %load configurations necessary for the script
    resultfolder = config.ResultFolder;
    type = config.type;
    
    %create storing folder for trajectory if not exist
    savefoldername = resultfolder+"/TwowayAnova_"+type+"/";
    if ~exist(savefoldername, 'dir')
       mkdir(savefoldername);
    end

    %processing the data into a long numeric vector 
    [YoungY, YoungGroupNames, YoungConditionNames]=ReGroupData(YoungData,'Young');
    [HealthyOldY, HealthyOldGroupNames, HealthyOldConditionNames]=ReGroupData(HealthyOldData,'HealthyOld');
    [MCIPosY, MCIPosGroupNames, MCIPosConditionNames]=ReGroupData(MCIPosData,'MCIPos');
    [MCINegY, MCINegGroupNames, MCINegConditionNames]=ReGroupData(MCINegData,'MCIPNeg');
    [MCIUnkY, MCIUnkGroupNames, MCIUnkConditionNames]=ReGroupData(MCIUnkData,'MCIPUnk');

    AllY = [MCIPosY,MCINegY,MCIUnkY,YoungY,HealthyOldY];
    AllGroupNames = [MCIPosGroupNames,MCINegGroupNames,MCIUnkGroupNames,YoungGroupNames,HealthyOldGroupNames];
    AllConditionNames = [MCIPosConditionNames,MCINegConditionNames,MCIUnkConditionNames,YoungConditionNames,HealthyOldConditionNames];

    %Do two-way anova with unbalanced design
    [p,anova_tab, stats]= anovan(AllY,{AllGroupNames,AllConditionNames},'model','interaction','varnames',{'Groups','Conditions'},'display','on');

    %Do multiple comparisons on main effect 1
    multicomp_tab1 = multcompare(stats,'Dimension',[1],'CType','bonferroni');
    title("Multiple comparisons with bonferroni correction of");
    saveas(gcf,savefoldername+"MultiCompME1.png");
    close(gcf);

    %Do multiple comparisons on main effect 2
    multicomp_tab2 = multcompare(stats,'Dimension',[2],'CType','bonferroni');
    title("Multiple comparisons with bonferroni correction");
    saveas(gcf,savefoldername+"MultiCompME2.png");
    close(gcf);

    %Do multiple comparisons on main effect 1&2
    multicomp_tab12 = multcompare(stats,'Dimension',[1,2],'CType','bonferroni');
    title("Multiple comparisons with bonferroni correction");
    saveas(gcf,savefoldername+"MultiCompME1ME2.png");
    close(gcf);    
end