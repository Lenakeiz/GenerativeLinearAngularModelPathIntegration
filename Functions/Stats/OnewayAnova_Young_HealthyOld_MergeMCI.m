%% One-way anova on merged MCI
function [anova_tab,multicomp_tab] = OnewayAnova_Young_HealthyOld_MergeMCI(YoungParams, HealthyOldParams, MCIParams, config)

    %load configurations necessary for the script
    resultfolder = config.ResultFolder;
    
    %create storing folder for trajectory if not exist
    savefoldername = resultfolder+"/OnewayAnova/";
    if ~exist(savefoldername, 'dir')
       mkdir(savefoldername);
    end

    ParamName = config.ParamName;
    param_nums = length(ParamName);    
    
    anova_tab = cell(0);
    multicomp_tab = cell(0);

    %for each parameter, calculate the anova and multiple test
    for param_idx=1:param_nums
        
        param_name = ParamName(param_idx);
        
        %processing the data into a long numeric vector 
        Params = [MCIParams(:,param_idx)',HealthyOldParams(:,param_idx)',YoungParams(:,param_idx)'];
        GroupNames = [string(repmat({'MCI'},1,size(MCIParams,1))),...
                      string(repmat({'HealthyOld'},1,size(HealthyOldParams,1))),...
                      string(repmat({'Young'},1,size(YoungParams,1)))];
    
        %Do one-way anova with unbalanced design
        [p,tb1, stats]= anova1(Params, GroupNames, 'display','on');
        anova_tab{param_idx} = tb1;
    
        %Do multiple comparisons on main effect 1
        result = multcompare(stats,'CType','bonferroni');
        multicomp_tab{param_idx} = result;
        title("Multiple comparisons with bonferroni correction of : "+param_name);
        saveas(gcf,savefoldername+"MultiComp"+param_name+".png");
        close(gcf);
    
    end
    
end
