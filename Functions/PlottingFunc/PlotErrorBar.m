function PlotErrorBar(YoungParamValues, HealthyControlsParamValues, MCIPosParamValues, MCINegParamValues, MCIUnkParamValues, config)
%   Error bar plot of estimated parameters across participants
%   YoungParamValues:
%   HealthyControlsParamValues: 
%   MCIPosParamValues:
%   MCINegParamValues:
%   MCIUnkParamValues:

%load configurations necessary for the script
TRIAL_FILTER = config.TrialFilter;
trialfolder = config.TrialFolder;

YoungParamValues_new = YoungParamValues;
YoungParamValues_new(:,end) = sqrt(1./YoungParamValues_new(:,end));

HealthyControlsParamValues_new = HealthyControlsParamValues;
HealthyControlsParamValues_new(:,end) = sqrt(1./HealthyControlsParamValues_new(:,end));

MCIPosParamValues_new = MCIPosParamValues;
MCIPosParamValues_new(:,end) = sqrt(1./MCIPosParamValues_new(:,end));

MCINegParamValues_new = MCINegParamValues;
MCINegParamValues_new(:,end) = sqrt(1./MCINegParamValues_new(:,end));

MCIUnkParamValues_new = MCIUnkParamValues;
MCIUnkParamValues_new(:,end) = sqrt(1./MCIUnkParamValues_new(:,end));


% Young
mean_Young = mean(YoungParamValues_new,1,'omitnan'); 
%std_Young = std(YoungParamValues,1);
sem_Young = std(YoungParamValues_new,1,'omitnan')./sqrt(size(YoungParamValues_new,1));
% Healthy
mean_Healthy = mean(HealthyControlsParamValues_new,1,'omitnan'); 
%std_Healthy = std(HealthyControlsParamValues,1);
sem_Healthy = std(HealthyControlsParamValues_new,1,'omitnan')./sqrt(size(HealthyControlsParamValues_new,1));
%MICPos
mean_MCIPos = mean(MCIPosParamValues_new,1,'omitnan');
%std_MCIPos = std(MCIPosParamValues,1);
sem_MCIPos = std(MCIPosParamValues_new,1,'omitnan')./sqrt(size(MCIPosParamValues_new,1));
%MCINeg
mean_MCINeg = mean(MCINegParamValues_new,1,'omitnan');
%std_MCINeg = std(MCINegParamValues,1);
sem_MCINeg = std(MCINegParamValues_new,1,'omitnan')./sqrt(size(MCINegParamValues_new,1));
%MCIUnk
mean_MCIUnk = mean(MCIUnkParamValues_new,1,'omitnan');
%std_MCIUnk = std(MCIUnkParamValues,1);
sem_MCIUnk = std(MCIUnkParamValues_new,1,'omitnan')./sqrt(size(MCIUnkParamValues_new,1));

%group data into parameter groups and each group has data from 5 kinds of
%people
mean_all = [mean_Young; mean_Healthy; mean_MCIPos; mean_MCINeg; mean_MCIUnk];
mean_all = mean_all'; %each row says each estimated parameters
%std_all = [std_Young; std_Healthy; std_MCIPos; std_MCINeg; std_MCIUnk];
%std_all = std_all'; %each row says each estimated parameters
sem_all = [sem_Young; sem_Healthy; sem_MCIPos; sem_MCINeg; sem_MCIUnk];
sem_all = sem_all';

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

hold on

%add y=1 for reference
yline(1);

hold on

%add y=pi/2 for reference
yline(pi/2, '--');

%add title
if TRIAL_FILTER==0
    title("All trials",'FontSize',20)
elseif TRIAL_FILTER==1
    title("No changes",'FontSize',20)
elseif TRIAL_FILTER==2
    title("No distal cues",'FontSize',20)
elseif TRIAL_FILTER==3
    title("No optic flow",'FontSize',20)
else
    error("Choose the correct trial flag!")
end

legend(b, {'Young' 'HealthyOld' 'MCIPos' 'MCINeg' 'MCIUnk'}, 'Location','northwest');
set(gca, 'XTickLabel', {'G_1' 'G_2' 'G_3' 'g_2' 'g_3' 'k_3' '\sigma' '\sqrt{1/\nu}'});    
exportgraphics(f,trialfolder+"errorbar.png",'Resolution',300);

hold off

end