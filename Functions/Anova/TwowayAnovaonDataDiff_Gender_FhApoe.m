%% Two-way anova on all Groups
% two groups:          AllApoePos / AllApoeNeg 
% three conditions:     no changed / no distal cue / no optical flow
function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnovaonDataDiff_Gender_FhApoe(PosPos_Diff,PosNeg_Diff, NegPos_Diff,NegNeg_Diff,config)

    %load configurations necessary for the script
    resultfolder = config.ResultFolder;
    
    %create storing folder for trajectory if not exist
    savefoldername = resultfolder+"/TwowayAnovaOnData_Gender_FhApoe/";
    if ~exist(savefoldername, 'dir')
       mkdir(savefoldername);
    end

    FHPos_ApoePos_Gender = config.FHPos_ApoePos_Gender;
    FHPos_ApoeNeg_Gender = config.FHPos_ApoeNeg_Gender;
    FHNeg_ApoePos_Gender = config.FHNeg_ApoePos_Gender;
    FHNeg_ApoeNeg_Gender = config.FHNeg_ApoeNeg_Gender;

    ParamName = config.ParamName;
    
    anova_tab = cell(0);
    multicomp_tab1 = cell(0);
    multicomp_tab2 = cell(0);
    multicomp_tab12 = cell(0);
    %for each parameter, calculate the anova and multiple test
        
    %processing the data into a long numeric vector 
    [PosPosY, PosPosGroupNames, PosPosGenderNames]=GroupAndRemoveNaN_2way_OnData(PosPos_Diff,'PosPos', FHPos_ApoePos_Gender);
    [PosNegY, PosNegGroupNames, PosNegGenderNames]=GroupAndRemoveNaN_2way_OnData(PosNeg_Diff,'PosNeg', FHPos_ApoeNeg_Gender);
    [NegPosY, NegPosGroupNames, NegPosGenderNames]=GroupAndRemoveNaN_2way_OnData(NegPos_Diff,'NegPos', FHNeg_ApoePos_Gender);
    [NegNegY, NegNegGroupNames, NegNegGenderNames]=GroupAndRemoveNaN_2way_OnData(NegNeg_Diff,'NegNeg', FHNeg_ApoeNeg_Gender);
    
    AllY = [PosPosY,PosNegY,NegPosY,NegNegY];
    AllGroupNames = [PosPosGroupNames,PosNegGroupNames,NegPosGroupNames,NegNegGroupNames];
    AllGenderNames = [PosPosGenderNames,PosNegGenderNames,NegPosGenderNames,NegNegGenderNames];

    %Do three-way anova with unbalanced design
    [p, tb1, stats]= anovan(AllY, {AllGroupNames,AllGenderNames}, 'model','interaction','varnames',{'Groups','Genders'},'display','on');
    anova_tab = tb1;

    %Do multiple comparisons on main effect 1
    result = multcompare(stats,'Dimension',[1],'CType','bonferroni');
    multicomp_tab1 = result;
    saveas(gcf,savefoldername+"MultiCompME1.png");
    %close(gcf);

    %Do multiple comparisons on main effect 2
    result = multcompare(stats,'Dimension',[2],'CType','bonferroni');
    multicomp_tab2 = result;
    saveas(gcf,savefoldername+"MultiCompME2.png");
    %close(gcf);

    %Do multiple comparisons on main effect 1&2
    result = multcompare(stats,'Dimension',[1,2],'CType','bonferroni');
    multicomp_tab12 = result;
    saveas(gcf,savefoldername+"MultiCompME1ME2.png");
    %close(gcf);    
    
end