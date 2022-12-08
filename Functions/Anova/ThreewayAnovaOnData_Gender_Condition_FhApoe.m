%% Three-way anova on the data
% two groups:           AllPos / AllNeg 
% three conditions:     no changed / no distal cue / no optical flow
% Genders:              1 male/ 2 female
function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab3, multicomp_tab12, multicomp_tab13, multicomp_tab23, multicomp_tab123] = ThreewayAnovaOnData_Gender_Condition_FhApoe(FHPos_ApoePos, FHPos_ApoeNeg, FHNeg_ApoePos, FHNeg_ApoeNeg, config)
    %FHPos_ApoePos_Params, FHPos_ApoeNeg_Param, FHNeg_ApoePos_Param, FHNeg_ApoeNeg_Param, config

    %load configurations necessary for the script
    resultfolder = config.ResultFolder;

    FHPos_ApoePos_Gender = config.FHPos_ApoePos_Gender;
    FHPos_ApoeNeg_Gender = config.FHPos_ApoeNeg_Gender;
    FHNeg_ApoePos_Gender = config.FHNeg_ApoePos_Gender;
    FHNeg_ApoeNeg_Gender = config.FHNeg_ApoeNeg_Gender;
    
    %create storing folder for trajectory if not exist
    savefoldername = resultfolder+"/ThreewayAnovaOnData_Gender_Condition_FhApoe/";
    if ~exist(savefoldername, 'dir')
       mkdir(savefoldername);
    end
    
    anova_tab = cell(0);
    multicomp_tab1 = cell(0);
    multicomp_tab2 = cell(0);
    multicomp_tab3 = cell(0);
    multicomp_tab12 = cell(0);
    multicomp_tab13 = cell(0);
    multicomp_tab23 = cell(0);
    multicomp_tab123 = cell(0);

        
        
    %processing the data into a long numeric vector 
    [PosPosY, PosPosGroupNames, PosPosConditionNames, PosPosGenderNames]=GroupAndRemoveNaN_3way_OnData(FHPos_ApoePos,'PosPos', FHPos_ApoePos_Gender);
    [PosNegY, PosNegGroupNames, PosNegConditionNames, PosNegGenderNames]=GroupAndRemoveNaN_3way_OnData(FHPos_ApoeNeg,'PosNeg', FHPos_ApoeNeg_Gender);
    [NegPosY, NegPosGroupNames, NegPosConditionNames, NegPosGenderNames]=GroupAndRemoveNaN_3way_OnData(FHNeg_ApoePos,'NegPos', FHNeg_ApoePos_Gender);
    [NegNegY, NegNegGroupNames, NegNegConditionNames, NegNegGenderNames]=GroupAndRemoveNaN_3way_OnData(FHNeg_ApoeNeg,'NegNeg', FHNeg_ApoeNeg_Gender);
    
    AllY = [PosPosY,PosNegY,NegPosY,NegNegY];
    AllGroupNames = [PosPosGroupNames,PosNegGroupNames,NegPosGroupNames,NegNegGroupNames];
    AllConditionNames = [PosPosConditionNames,PosNegConditionNames,NegPosConditionNames,NegNegConditionNames];
    AllGenderNames = [PosPosGenderNames, PosNegGenderNames,NegPosGenderNames,NegNegGenderNames];

    %Do three-way anova with unbalanced design
    [p, tb1, stats]= anovan(AllY, {AllGroupNames,AllConditionNames, AllGenderNames}, 'model','full','varnames',{'Groups','Conditions','Genders'},'display','on');
    anova_tab = tb1;

    %Do multiple comparisons on main effect 1
    figure()
    result = multcompare(stats,'Dimension',[1],'CType','bonferroni');
    multicomp_tab1 = result;
    saveas(gcf,savefoldername+"ME1.png");
    close(gcf);

    %Do multiple comparisons on main effect 2
    figure()
    result = multcompare(stats,'Dimension',[2],'CType','bonferroni');
    multicomp_tab2 = result;
    saveas(gcf,savefoldername+"ME2.png");
    close(gcf);

    %Do multiple comparisons on main effect 3
    figure()
    result = multcompare(stats,'Dimension',[3],'CType','bonferroni');
    multicomp_tab3 = result;
    saveas(gcf,savefoldername+"ME3.png");
    close(gcf);

    %Do multiple comparisons on main effect 1&2
    figure()
    result = multcompare(stats,'Dimension',[1 2],'CType','bonferroni');
    multicomp_tab12 = result;
    saveas(gcf,savefoldername+"ME1ME2.png");
    close(gcf);    

    %Do multiple comparisons on main effect 1&3
    figure()
    result = multcompare(stats,'Dimension',[1 3],'CType','bonferroni');
    multicomp_tab13 = result;
    saveas(gcf,savefoldername+"ME1ME3.png");
    close(gcf); 

    %Do multiple comparisons on main effect 2&3
    figure()
    result = multcompare(stats,'Dimension',[2 3],'CType','bonferroni');
    multicomp_tab23 = result;
    saveas(gcf,savefoldername+"ME2ME3.png");
    close(gcf); 

    %Do multiple comparisons on main effect 1&2&3
    figure()
    result = multcompare(stats,'Dimension',[1 2 3],'CType','bonferroni');
    multicomp_tab123 = result;
    saveas(gcf,savefoldername+"ME1ME2ME3.png");
    close(gcf);  
    
end