%% Two-way anova on all Groups
% two groups:          AllApoePos / AllApoeNeg 
% three conditions:     no changed / no distal cue / no optical flow
function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnovaonParameterDiff_Gender_FhApoe(PosPos_Param_Diff,PosNeg_Param_Diff, NegPos_Param_Diff,NegNeg_Param_Diff,config)

    %load configurations necessary for the script
    resultfolder = config.ResultFolder;
    
    %create storing folder for trajectory if not exist
    savefoldername = resultfolder+"/TwowayAnovaOnParameterDiff_Gender_FhApoe/";
    if ~exist(savefoldername, 'dir')
       mkdir(savefoldername);
    end

    FHPos_ApoePos_Gender = config.FHPos_ApoePos_Gender;
    FHPos_ApoeNeg_Gender = config.FHPos_ApoeNeg_Gender;
    FHNeg_ApoePos_Gender = config.FHNeg_ApoePos_Gender;
    FHNeg_ApoeNeg_Gender = config.FHNeg_ApoeNeg_Gender;

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
        [PosPosY, PosPosGroupNames, PosPosGenderNames]=GroupAndRemoveNaN_ParameterDiff(PosPos_Param_Diff,param_idx,'PosPos', FHPos_ApoePos_Gender);
        [PosNegY, PosNegGroupNames, PosNegGenderNames]=GroupAndRemoveNaN_ParameterDiff(PosNeg_Param_Diff,param_idx,'PosNeg', FHPos_ApoeNeg_Gender);
        [NegPosY, NegPosGroupNames, NegPosGenderNames]=GroupAndRemoveNaN_ParameterDiff(NegPos_Param_Diff,param_idx,'NegPos', FHNeg_ApoePos_Gender);
        [NegNegY, NegNegGroupNames, NegNegGenderNames]=GroupAndRemoveNaN_ParameterDiff(NegNeg_Param_Diff,param_idx,'NegNeg', FHNeg_ApoeNeg_Gender);
        
        AllY = [PosPosY,PosNegY,NegPosY,NegNegY];
        AllGroupNames = [PosPosGroupNames,PosNegGroupNames,NegPosGroupNames,NegNegGroupNames];
        AllGenderNames = [PosPosGenderNames,PosNegGenderNames,NegPosGenderNames,NegNegGenderNames];
    
        %Do three-way anova with unbalanced design
        [p, tb1, stats]= anovan(AllY, {AllGroupNames,AllGenderNames}, 'model','interaction','varnames',{'Groups','Genders'},'display','on');
        anova_tab{param_idx} = tb1;
    
        %Do multiple comparisons on main effect 1
        result = multcompare(stats,'Dimension',[1],'CType','bonferroni');
        multicomp_tab1{param_idx} = result;
        title("Parameter: "+param_name);
        saveas(gcf,savefoldername+"MultiCompME1_"+param_name+".png");
        %close(gcf);
    
        %Do multiple comparisons on main effect 2
        result = multcompare(stats,'Dimension',[2],'CType','bonferroni');
        multicomp_tab2{param_idx} = result;
        title("Parameter: "+param_name);
        saveas(gcf,savefoldername+"MultiCompME2_"+param_name+".png");
        %close(gcf);
    
        %Do multiple comparisons on main effect 1&2
        result = multcompare(stats,'Dimension',[1,2],'CType','bonferroni');
        multicomp_tab12{param_idx} = result;
        title("Parameter: "+param_name);
        saveas(gcf,savefoldername+"MultiCompME1ME2_"+param_name+".png");
        %close(gcf);    
    end
    
end