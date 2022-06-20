%% Two-way anova on all Groups
% five groups:          Young / HealthyOld / MCIPos / MCINeg / MCIUnk
% three conditions:     no changed / no distal cue / no optical flow
function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab12] = TwowayAnova_MCIPosMCINeg(AllMCIPos, AllMCINeg, config)

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
        [MCIPosY, MCIPosGroupNames, MCIPosConditionNames]=GroupAndRemoveNaN(AllMCIPos,param_idx,'MCIPos');
        [MCINegY, MCINegGroupNames, MCINegConditionNames]=GroupAndRemoveNaN(AllMCINeg,param_idx,'MCIPNeg');
    
        AllY = [MCIPosY,MCINegY];
        AllGroupNames = [MCIPosGroupNames,MCINegGroupNames];
        AllConditionNames = [MCIPosConditionNames,MCINegConditionNames];
    
        %Do two-way anova with unbalanced design
        [p,tb1, stats]= anovan(AllY,{AllGroupNames,AllConditionNames},'model','interaction','varnames',{'Groups','Conditions'},'display','on');
        anova_tab{param_idx} = tb1;
    
        %Do multiple comparisons on main effect 1
        result = multcompare(stats,'Dimension',[1],'CType','bonferroni');
        multicomp_tab1{param_idx} = result;
        title("parameter: "+param_name);
        saveas(gcf,savefoldername+"MultiCompME1_"+param_name+".png");
        close(gcf);
    
        %Do multiple comparisons on main effect 2
        result = multcompare(stats,'Dimension',[2],'CType','bonferroni');
        multicomp_tab2{param_idx} = result;
        title("parameter: "+param_name);
        saveas(gcf,savefoldername+"MultiCompME2_"+param_name+".png");
        close(gcf);
    
        %Do multiple comparisons on main effect 1&2
        result = multcompare(stats,'Dimension',[1,2],'CType','bonferroni');
        multicomp_tab12{param_idx} = result;
        title("parameter: "+param_name);
        saveas(gcf,savefoldername+"MultiCompME1ME2_"+param_name+".png");
        close(gcf);    
    end
    
end
