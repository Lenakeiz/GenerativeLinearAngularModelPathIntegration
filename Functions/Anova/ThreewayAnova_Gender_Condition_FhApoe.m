%% Three-way anova on the data
% two groups:           AllPos / AllNeg 
% three conditions:     no changed / no distal cue / no optical flow
% Genders:              1 male/ 2 female
function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab3, multicomp_tab12, multicomp_tab13, multicomp_tab23, multicomp_tab123] = ThreewayAnova_Gender_Condition_FhApoe(FHPos_ApoePos_Params, FHPos_ApoeNeg_Param, FHNeg_ApoePos_Param, FHNeg_ApoeNeg_Param, config)
    %FHPos_ApoePos_Params, FHPos_ApoeNeg_Param, FHNeg_ApoePos_Param, FHNeg_ApoeNeg_Param, config

    %load configurations necessary for the script
    resultfolder = config.ResultFolder;

    FHPos_ApoePos_Gender = config.FHPos_ApoePos_Gender;
    FHPos_ApoeNeg_Gender = config.FHPos_ApoeNeg_Gender;
    FHNeg_ApoePos_Gender = config.FHNeg_ApoePos_Gender;
    FHNeg_ApoeNeg_Gender = config.FHNeg_ApoeNeg_Gender;
    
    %create storing folder for trajectory if not exist
    savefoldername = resultfolder+"/ThreewayAnova_Gender_Condition_FhApoe/";
    if ~exist(savefoldername, 'dir')
       mkdir(savefoldername);
    end
    
    ParamName = config.ParamName;
    param_nums = length(ParamName);
    
    anova_tab = cell(0);
    multicomp_tab1 = cell(0);
    multicomp_tab2 = cell(0);
    multicomp_tab3 = cell(0);
    multicomp_tab12 = cell(0);
    multicomp_tab13 = cell(0);
    multicomp_tab23 = cell(0);
    multicomp_tab123 = cell(0);

    %for each parameter, calculate the anova and multiple test
    for param_idx=1:param_nums
        
        param_name = ParamName(param_idx);
        
        %processing the data into a long numeric vector 
        [PosPosY, PosPosGroupNames, PosPosConditionNames, PosPosGenderNames]=GroupAndRemoveNaN_3way(FHPos_ApoePos_Params,param_idx,'PosPos', FHPos_ApoePos_Gender);
        [PosNegY, PosNegGroupNames, PosNegConditionNames, PosNegGenderNames]=GroupAndRemoveNaN_3way(FHPos_ApoeNeg_Param,param_idx,'PosNeg', FHPos_ApoeNeg_Gender);
        [NegPosY, NegPosGroupNames, NegPosConditionNames, NegPosGenderNames]=GroupAndRemoveNaN_3way(FHNeg_ApoePos_Param,param_idx,'NegPos', FHNeg_ApoePos_Gender);
        [NegNegY, NegNegGroupNames, NegNegConditionNames, NegNegGenderNames]=GroupAndRemoveNaN_3way(FHNeg_ApoeNeg_Param,param_idx,'NegNeg', FHNeg_ApoeNeg_Gender);
        
        AllY = [PosPosY,PosNegY,NegPosY,NegNegY];
        AllGroupNames = [PosPosGroupNames,PosNegGroupNames,NegPosGroupNames,NegNegGroupNames];
        AllConditionNames = [PosPosConditionNames,PosNegConditionNames,NegPosConditionNames,NegNegConditionNames];
        AllGenderNames = [PosPosGenderNames, PosNegGenderNames,NegPosGenderNames,NegNegGenderNames];
    
        %Do three-way anova with unbalanced design
        [p, tb1, stats]= anovan(AllY, {AllGroupNames,AllConditionNames, AllGenderNames}, 'model','full','varnames',{'Groups','Conditions','Genders'},'display','on');
        anova_tab{param_idx} = tb1;
    
        %Do multiple comparisons on main effect 1
        figure()
        result = multcompare(stats,'Dimension',[1],'CType','bonferroni');
        multicomp_tab1{param_idx} = result;
        title("Parameter: "+param_name);
        saveas(gcf,savefoldername+param_name+"_ME1.png");
        close(gcf);
    
        %Do multiple comparisons on main effect 2
        figure()
        result = multcompare(stats,'Dimension',[2],'CType','bonferroni');
        multicomp_tab2{param_idx} = result;
        title("Parameter: "+param_name);
        saveas(gcf,savefoldername+param_name+"_ME2.png");
        close(gcf);
 
        %Do multiple comparisons on main effect 3
        figure()
        result = multcompare(stats,'Dimension',[3],'CType','bonferroni');
        multicomp_tab3{param_idx} = result;
        title("Parameter: "+param_name);
        saveas(gcf,savefoldername+param_name+"_ME3.png");
        close(gcf);

        %Do multiple comparisons on main effect 1&2
        figure()
        result = multcompare(stats,'Dimension',[1 2],'CType','bonferroni');
        multicomp_tab12{param_idx} = result;
        title("Parameter: "+param_name);
        saveas(gcf,savefoldername+param_name+"_ME1ME2.png");
        close(gcf);    

        %Do multiple comparisons on main effect 1&3
        figure()
        result = multcompare(stats,'Dimension',[1 3],'CType','bonferroni');
        multicomp_tab13{param_idx} = result;
        title("Parameter: "+param_name);
        saveas(gcf,savefoldername+param_name+"_ME1ME3.png");
        close(gcf); 

        %Do multiple comparisons on main effect 2&3
        figure()
        result = multcompare(stats,'Dimension',[2 3],'CType','bonferroni');
        multicomp_tab23{param_idx} = result;
        title("Parameter: "+param_name);
        saveas(gcf,savefoldername+param_name+"_ME2ME3.png");
        close(gcf); 

        %Do multiple comparisons on main effect 1&2&3
        figure()
        result = multcompare(stats,'Dimension',[1 2 3],'CType','bonferroni');
        multicomp_tab123{param_idx} = result;
        title("Parameter: "+param_name);
        saveas(gcf,savefoldername+param_name+"_ME1ME2ME3.png");
        close(gcf);  
    end
    
end