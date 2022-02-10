function TwowayAnovaOnParamsMergeMCI(AllYoung, AllHealthyOld, AllMCI, config)
%   Error bar plot of estimated parameters across participants
%   AllYoung
%   AllHealthyOld
%   AllMCI

%load configurations necessary for the script
resultfolder = config.ResultFolder;

%create storing folder for trajectory if not exist
savefoldername = resultfolder+"TwowayAnovaOnParamsMergeMCI/";
if ~exist(savefoldername, 'dir')
   mkdir(savefoldername);
end

param_names = ["gamma", "bG3", "g2", "g3", 'b', "sigma", "nu"];
param_nums = length(param_names);

anova_tab = cell(0);
multicomp_tab = cell(0);
%for each parameter, calculate the anova and multiple test
for param_idx=1:param_nums
    
    param_name = param_names(param_idx);
    
    %processing the data into a long numeric vector 
    [MCIY, MCIGroupNames, MCIConditionNames]=GroupAndRemoveNaN(AllMCI,param_idx,'MCI');
    [HealthyOldY, HealthyOldGroupNames, HealthyOldConditionNames]=GroupAndRemoveNaN(AllHealthyOld,param_idx,'HealthyOld');
    [YoungY, YoungGroupNames, YoungConditionNames]=GroupAndRemoveNaN(AllYoung,param_idx,'Young');
    
    AllY = [MCIY,HealthyOldY,YoungY];
    AllGroupNames = [MCIGroupNames,HealthyOldGroupNames,YoungGroupNames];
    AllConditionNames = [MCIConditionNames,HealthyOldConditionNames,YoungConditionNames];

    %Do two-way anova with unbalanced design
    [p,tb1, stats]= anovan(AllY,{AllGroupNames,AllConditionNames},'model','interaction','varnames',{'Groups','Conditions'},'display','on');
    anova_tab{param_idx} = tb1;
%     title("Standard ANOVA table of parameter: "+param_name);
%     saveas(gcf,savefoldername+"TwowayAnova_"+param_name+".png");
%     close(gcf)

    %Do multiple comparisons on main effect 1
    result = multcompare(stats,'Dimension',[1],'CType','bonferroni');
    multicomp_tab{param_idx} = result;
    title("Multiple comparisons with bonferroni correction of : "+param_name);
    saveas(gcf,savefoldername+"MultiCompME1_"+param_name+".png");
    close(gcf);

    %Do multiple comparisons on main effect 2
    result = multcompare(stats,'Dimension',[2],'CType','bonferroni');
    multicomp_tab{param_idx} = result;
    title("Multiple comparisons with bonferroni correction of parameter: "+param_name);
    saveas(gcf,savefoldername+"MultiCompME2_"+param_name+".png");
    close(gcf);

    %Do multiple comparisons on main effect 1&2
    result = multcompare(stats,'Dimension',[1,2],'CType','bonferroni');
    multicomp_tab{param_idx} = result;
    title("Multiple comparisons with bonferroni correction of parameter: "+param_name);
    saveas(gcf,savefoldername+"MultiCompME1ME2_"+param_name+".png");
    close(gcf);    
end

save(savefoldername+"anova_info.mat", 'anova_tab');
save(savefoldername+"multicomp_info.mat", 'multicomp_tab');

end