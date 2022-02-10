%% model comparision
% Cleaning variables
clearvars; clear all; close all;

%load the BIC for base model
BIC_base_add_G3 = load("C:/Users/Zilong/Desktop/path integration model/model_estimated_x3/Output/BIC_base_add_G3.mat");
mean_base_add_G3 = zeros(1,5); sem_base_add_G3 = zeros(1,5);
for i=1:5 %for five different groups
    mean_base_add_G3(i) = mean(BIC_base_add_G3.BIC_ALL{i});
    sem_base_add_G3(i) = std(BIC_base_add_G3.BIC_ALL{i})./sqrt(size(BIC_base_add_G3.BIC_ALL{i},2));
end

%load the BIC for base model
BIC_base = load("C:/Users/Zilong/Desktop/path integration model/model_estimated_x3/Output/BIC_base.mat");
mean_base = zeros(1,5); sem_base = zeros(1,5);
for i=1:5 %for five different groups
    mean_base(i) = mean(BIC_base.BIC_ALL{i});
    sem_base(i) = std(BIC_base.BIC_ALL{i})./sqrt(size(BIC_base.BIC_ALL{i},2));
end

%load the BIC for g3=1 model
BIC_set_k3_0 = load("C:/Users/Zilong/Desktop/path integration model/model_estimated_x3/Output/BIC_set_k3_0.mat");
mean_set_k3_0= zeros(1,5); sem_set_k3_0 = zeros(1,5);
for i=1:5 %for five different groups
    mean_set_k3_0(i) = mean(BIC_set_k3_0.BIC_ALL{i});
    sem_set_k3_0(i) = std(BIC_set_k3_0.BIC_ALL{i})./sqrt(size(BIC_set_k3_0.BIC_ALL{i},2));
end

%load the BIC for g3=1 model
BIC_set_g3_1_k3_0 = load("C:/Users/Zilong/Desktop/path integration model/model_estimated_x3/Output/BIC_set_g3_1_k3_0.mat");
mean_set_g3_1_k3_0= zeros(1,5); sem_set_g3_1_k3_0 = zeros(1,5);
for i=1:5 %for five different groups
    mean_set_g3_1_k3_0(i) = mean(BIC_set_g3_1_k3_0.BIC_ALL{i});
    sem_set_g3_1_k3_0(i) = std(BIC_set_g3_1_k3_0.BIC_ALL{i})./sqrt(size(BIC_set_g3_1_k3_0.BIC_ALL{i},2));
end

%load the BIC for g2=1 model
BIC_set_g2_1 = load("C:/Users/Zilong/Desktop/path integration model/model_estimated_x3/Output/BIC_set_g2_1.mat");
mean_set_g2_1 = zeros(1,5); sem_set_g2_1 = zeros(1,5);
for i=1:5 %for five different groups
    mean_set_g2_1(i) = mean(BIC_set_g2_1.BIC_ALL{i});
    sem_set_g2_1(i) = std(BIC_set_g2_1.BIC_ALL{i})./sqrt(size(BIC_set_g2_1.BIC_ALL{i},2));
end

%load the BIC for G2=1 model
BIC_set_bG2_1 = load("C:/Users/Zilong/Desktop/path integration model/model_estimated_x3/Output/BIC_set_bG2_1.mat");
mean_set_bG2_1 = zeros(1,5); sem_set_bG2_1 = zeros(1,5);
for i=1:5 %for five different groups
    mean_set_bG2_1(i) = mean(BIC_set_bG2_1.BIC_ALL{i});
    sem_set_bG2_1(i) = std(BIC_set_bG2_1.BIC_ALL{i})./sqrt(size(BIC_set_bG2_1.BIC_ALL{i},2));
end

%load the BIC for G1=1 model
BIC_set_bG1_1 = load("C:/Users/Zilong/Desktop/path integration model/model_estimated_x3/Output/BIC_set_bG1_1.mat");
mean_set_bG1_1 = zeros(1,5); sem_set_bG1_1 = zeros(1,5);
for i=1:5 %for five different groups
    mean_set_bG1_1(i) = mean(BIC_set_bG1_1.BIC_ALL{i});
    sem_set_bG1_1(i) = std(BIC_set_bG1_1.BIC_ALL{i})./sqrt(size(BIC_set_bG1_1.BIC_ALL{i},2));
end

%load the BIC for G1=G2 model
BIC_set_same_bG = load("C:/Users/Zilong/Desktop/path integration model/model_estimated_x3/Output/BIC_set_same_bG.mat");
mean_set_same_bG = zeros(1,5); sem_set_BIC_set_same_bG = zeros(1,5);
for i=1:5 %for five different groups
    mean_set_same_bG(i) = mean(BIC_set_same_bG.BIC_ALL{i});
    sem_set_BIC_set_same_bG(i) = std(BIC_set_same_bG.BIC_ALL{i})./sqrt(size(BIC_set_same_bG.BIC_ALL{i},2));
end

%load the BIC for g2=g3 model
BIC_set_same_sg = load("C:/Users/Zilong/Desktop/path integration model/model_estimated_x3/Output/BIC_set_same_sg.mat");
mean_set_same_sg = zeros(1,5); sem_set_BIC_set_same_sg = zeros(1,5);
for i=1:5 %for five different groups
    mean_set_same_sg(i) = mean(BIC_set_same_sg.BIC_ALL{i});
    sem_set_BIC_set_same_sg(i) = std(BIC_set_same_sg.BIC_ALL{i})./sqrt(size(BIC_set_same_sg.BIC_ALL{i},2));
end

%load the BIC for g2=g3 model
BIC_set_same_Gg = load("C:/Users/Zilong/Desktop/path integration model/model_estimated_x3/Output/BIC_set_same_Gg.mat");
mean_set_same_Gg = zeros(1,5); sem_set_BIC_set_same_Gg = zeros(1,5);
for i=1:5 %for five different groups
    mean_set_same_Gg(i) = mean(BIC_set_same_Gg.BIC_ALL{i});
    sem_set_BIC_set_same_Gg(i) = std(BIC_set_same_Gg.BIC_ALL{i})./sqrt(size(BIC_set_same_Gg.BIC_ALL{i},2));
end

mean_all = [mean_base_add_G3; 
            mean_base; 
            mean_set_k3_0; mean_set_g3_1_k3_0; mean_set_g2_1; mean_set_bG2_1; mean_set_bG1_1; 
            mean_set_same_bG; mean_set_same_sg; mean_set_same_Gg];
mean_all = mean_all'; %each row says each group each col says each model
sem_all = [sem_base_add_G3; 
           sem_base; 
           sem_set_k3_0; sem_set_g3_1_k3_0; sem_set_g2_1; sem_set_bG2_1; sem_set_bG1_1; 
           sem_set_BIC_set_same_bG; sem_set_BIC_set_same_sg; sem_set_BIC_set_same_Gg];
sem_all = sem_all'; %each row says each group each col says each model

f = figure('visible','off','Position', [100 100 1000 800]);

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

legend(b, {'Base-add-G3', ...
           'Base' ...
           'g_3=1' 'g_3=1,k_3=0' 'g_2=1' 'G_2=1' 'G_1=1' ...
           'G_1=G_2' 'g_2=g_3' 'G_1=G_2,g_2=g_3'}, ...
           'Location','northoutside','FontSize',10, 'NumColumns', 5);

ylabel('Bayesian Inference Criterion (BIC)', 'FontSize',20);
%ylabel('Negative Loglikelihood','FontSize',20);

ylim([4,7]);

set(gca, 'XTickLabel', {'Young' 'Healthy' 'MCIPos' 'MCINeg' 'MCIUnk'});

exportgraphics(f,"C:/Users/Zilong/Desktop/path integration model/model_estimated_x3/Output/BIC.png",'Resolution',300);


