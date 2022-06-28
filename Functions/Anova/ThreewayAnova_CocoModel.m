%% Three-way anova on the data
% two groups:           AllPos / AllNeg 
% three conditions:     no changed / no distal cue / no optical flow
% Genders:              1 male/ 2 female
function [anova_tab,multicomp_tab1,multicomp_tab2, multicomp_tab3, multicomp_tab12, multicomp_tab13, multicomp_tab23, multicomp_tab123] = ThreewayAnova_CocoModel(AllPos, AllNeg, config)

    %load configurations necessary for the script
    resultfolder = config.ResultFolder;

    PosGender = config.Gender_Pos;
    NegGender = config.Gender_Neg;
    
    %create storing folder for trajectory if not exist
    savefoldername = resultfolder+"/ThreewayAnova/";
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
        [PosY, PosGroupNames, PosConditionNames, PosGenderNames]=GroupAndRemoveNaN_3way(AllPos,param_idx,'Pos', PosGender);
        [NegY, NegGroupNames, NegConditionNames, NegGenderNames]=GroupAndRemoveNaN_3way(AllNeg,param_idx,'Neg', NegGender);
    
        AllY = [PosY,NegY];
        AllGroupNames = [PosGroupNames,NegGroupNames];
        AllConditionNames = [PosConditionNames,NegConditionNames];
        AllGenderNames = [PosGenderNames, NegGenderNames];
    
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