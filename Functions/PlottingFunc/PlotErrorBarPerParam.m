function PlotErrorBarPerParam(AllYoung, AllHealthyOld, AllMCIPos, AllMCINeg, AllMCIUnk, config)
%   Error bar plot of estimated parameters across participants
%   YoungParamValues:
%   HealthyControlsParamValues: 
%   MCIPosParamValues:
%   MCINegParamValues:
%   MCIUnkParamValues:

resultfolder = config.ResultFolder;

%create storing folder for trajectory if not exist
savefoldername = resultfolder+"ErrorBar/";
if ~exist(savefoldername, 'dir')
   mkdir(savefoldername);
end

title_table = ["\gamma", "bG_3", "g_2", "g_3", "k_3", "\sigma", "\nu"];

numConds = size(AllYoung,2);
numParams = size(AllYoung{1},2);

mean_all = zeros(numConds, 5); %5 is the group number
sem_all = zeros(numConds, 5);

for paramidx = 1:numParams

    for trial_filter = 1:numConds
        YoungParam = AllYoung{trial_filter}(:,paramidx);
        [young_m, young_s] = getMeanSem(YoungParam);
        mean_all(trial_filter,1) = young_m; sem_all(trial_filter,1)=young_s;

        HealthyOldParam = AllHealthyOld{trial_filter}(:,paramidx);
        [healthy_m, healthy_s] = getMeanSem(HealthyOldParam);
        mean_all(trial_filter,2) = healthy_m; sem_all(trial_filter,2)=healthy_s;

        MCIPosParam = AllMCIPos{trial_filter}(:,paramidx);
        [mcipos_m, mcipos_s] = getMeanSem(MCIPosParam);
        mean_all(trial_filter,3) = mcipos_m; sem_all(trial_filter,3)=mcipos_s;

        MCINegParam = AllMCINeg{trial_filter}(:,paramidx);
        [mcineg_m, mcineg_s] = getMeanSem(MCINegParam);
        mean_all(trial_filter,4) = mcineg_m; sem_all(trial_filter,4)=mcineg_s;

        MCIUnkParam = AllMCIUnk{trial_filter}(:,paramidx);
        [mciunk_m, mciunk_s] = getMeanSem(MCIUnkParam);
        mean_all(trial_filter,5) = mciunk_m; sem_all(trial_filter,5)=mciunk_s;        
    end

    %create a figure for each parameter for plotting error bar
    f = figure('visible','off','Position', [100 100 1000 500]);
    ylabel('Estimated value of parameters');
    
    b = bar(mean_all, 'grouped');
    hold on

    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(mean_all);
    % Get the x coordinate of the bars
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end
    % Plot the errorbars
    errorbar(x',mean_all,sem_all,'k','linestyle','none');

    if title_table(paramidx)=="k_3" 
        hold on
        %add y=pi/2 for angular reference
        yline(pi/2, '--');
        yline(pi);
    else
        %add reference lines
        hold on
        %add y=1 for reference
        yline(1);
    end

    %add title
    title(title_table(paramidx),'FontSize',20)

    legend(b, {'Young' 'HealthyOld' 'MCIPos' 'MCINeg' 'MCIUnk'}, 'Location','northwest', 'NumColumns',5);
    set(gca, 'XTickLabel', {'No change' 'No distal cue' 'No optical flow'});    
    exportgraphics(f,savefoldername+title_table(paramidx)+".png",'Resolution',300);
    
    hold off

end

function [m, s] = getMeanSem(Param)
    m = mean(Param);
    s = std(Param)./sqrt(length(Param));
end

end
