function OnewayAnovaOnParams(YoungParamValues, HealthyControlsParamValues, MCIPosParamValues, MCINegParamValues, MCIUnkParamValues, config)
%   Error bar plot of estimated parameters across participants
%   YoungParamValues:
%   HealthyControlsParamValues: 
%   MCIPosParamValues:
%   MCINegParamValues:
%   MCIUnkParamValues:

%load configurations necessary for the script
trialfolder = config.TrialFolder;

%create storing folder for trajectory if not exist
savefoldername = trialfolder+"OnewayAnova/";
if ~exist(savefoldername, 'dir')
   mkdir(savefoldername);
end

param_names = ["bG1", "bG2", "bG3", "g2", "g3", 'k3', "sigma", "nu"];
param_nums = size(YoungParamValues,2);
anova_tab = cell(0);
multicomp_tab = cell(0);
%for each parameter, calculate the anova and multiple test
for param_idx=1:param_nums
    
    param_name = param_names(param_idx);

    param_values = [MCIPosParamValues(:,param_idx);
                    MCINegParamValues(:,param_idx);
                    MCIUnkParamValues(:,param_idx);
                    YoungParamValues(:,param_idx);
                    HealthyControlsParamValues(:,param_idx)];
    param_values = param_values';
    
    group_names = [string(repmat({'MCIPos'},size(MCIPosParamValues,1),1));
                   string(repmat({'MCINeg'},size(MCINegParamValues,1),1));
                   string(repmat({'MCIUnk'},size(MCIUnkParamValues,1),1));
                   string(repmat({'Young'},size(YoungParamValues,1),1));
                   string(repmat({'OldHealthy'},size(HealthyControlsParamValues,1),1))];
    group_names = group_names';
    
    %Do one-way anova on the estimated parameter
    [p, tb, stats] = anova1(param_values,group_names);
    anova_tab{param_idx} = tb;
    title("Box plot of parameter: "+param_name);
    saveas(gcf,savefoldername+"Anova_"+param_name+".png");
    close(gcf);

    %Do multiple comparisons
    %result = multcompare(stats);
    result = multcompare(stats,'CType','bonferroni');
    multicomp_tab{param_idx} = result;
    title("multiple comparisons with bonferroni correction of parameter: "+param_name);
    saveas(gcf,savefoldername+"MultiComp_"+param_name+".png");
    close(gcf);
end

save(savefoldername+"anova_info.mat", 'anova_tab');
save(savefoldername+"multicomp_info.mat", 'multicomp_tab');

end