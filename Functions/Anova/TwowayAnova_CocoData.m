%% Two-way anova on all Groups
% two groups:          AllApoePos / AllApoeNeg 
% three conditions:     no changed / no distal cue / no optical flow
function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnova_CocoData(AllPos, AllNeg, config)

    %load configurations necessary for the script
    resultfolder = config.ResultFolder;
    
    %create storing folder for trajectory if not exist
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
        [PosY, PosGroupNames, PosConditionNames]=GroupAndRemoveNaN(AllPos,param_idx,'Pos');
        [NegY, NegGroupNames, NegConditionNames]=GroupAndRemoveNaN(AllNeg,param_idx,'Neg');
    
        AllY = [PosY,NegY];
        AllGroupNames = [PosGroupNames,NegGroupNames];
        AllConditionNames = [PosConditionNames,NegConditionNames];
    
        %Do three-way anova with unbalanced design
        [p, tb1, stats]= anovan(AllY, {AllGroupNames,AllConditionNames}, 'model','interaction','varnames',{'Groups','Conditions'},'display','on');
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