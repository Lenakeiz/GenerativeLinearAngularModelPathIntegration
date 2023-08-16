%% Script to check the age and education difference in MCI+ and MCI-
% Zilong Ji, UCL, 2022 zilong.ji@ucl.ac.uk

% Preparing the data
GLAMPI_PrepareBaseConfig

% Preprocessing the data
GLAMPI_PreprocessData

%load('Data/HowettBrain2019_Dataset.mat');

%%
ColorPattern; 

MCINeg_ID = {'A093863_13', 'A093863_35', 'A093863_42', 'A093863_47', 'A093863_53', 'A093863_55', 'A093863_56', 'A093863_66', 'A093863_69', 'A093863_75', 'A093863_81', 'A093863_08', 'A093863_04', 'A093863_50'};

%IDs not participate in ANOVA analysis are excluded, that is the last two,
%so 9 IDs left
MCIPos_ID = {'A093863_21', 'A093863_27', 'A093863_37',	'A093863_39', 'A093863_41',	'A093863_62', 'A093863_68',	'A093863_59', 'A093863_70', 'A093863_40', 'A093863_49'};

HO_ID = {'A093863_11',	'A093863_12',	'A093863_15',	'A093863_16',	'A093863_18',	'A093863_20',	'A093863_26',	'A093863_29',	'A093863_31',	'A093863_34',	'A093863_36',	'A093863_38',	'A093863_43',	'A093863_44',	'A093863_51',	'A093863_52',	'A093863_54',	'A093863_57',	'A093863_60',	'A093863_61',	'A093863_06',	'A093863_71',	'A093863_73',	'A093863_74', 'A093863_77',	'A093863_78',	'A093863_80',	'A093863_83',	'A093863_84',	'A093863_88', 'A093863_09', 'A093863_22', 'A093863_23'};

%load csv
PsychTab = readtable('C:\Users\zji\Desktop\DavidData\HowettBrain2019_NeuroPsych.csv');

IDs = PsychTab.ID;
Ages = PsychTab.Age;
YrsinEds = PsychTab.YrsinEd;

%mcineg
mcineg_age = zeros(1,size(MCINeg_ID,2));
mcineg_ed = zeros(1,size(MCINeg_ID,2));
for i=1:size(MCINeg_ID,2)
    rowIndex = find(strcmp(IDs, MCINeg_ID{i}));
    mcineg_age(i) = Ages(rowIndex);
    mcineg_ed(i) = YrsinEds(rowIndex);
end

%mcipos
mcipos_age = zeros(1,size(MCIPos_ID,2));
mcipos_ed = zeros(1,size(MCIPos_ID,2));
for i=1:size(MCIPos_ID,2)
    rowIndex = find(strcmp(IDs, MCIPos_ID{i}));
    mcipos_age(i) = Ages(rowIndex);
    mcipos_ed(i) = YrsinEds(rowIndex);
end

%healthy old
ho_age = zeros(1,size(HO_ID,2));
ho_ed = zeros(1,size(HO_ID,2));
for i=1:size(HO_ID,2)
    rowIndex = find(strcmp(IDs, HO_ID{i}));
    ho_age(i) = Ages(rowIndex);
    ho_ed(i) = YrsinEds(rowIndex);
end

%make dir
ResultFolder = pwd + "/Output/FigS7";

if ~exist(ResultFolder, 'dir')
   mkdir(ResultFolder);
end


%% compare the age between MCI+ and MCI-
f = figure('Position', [100 100 600 300]);

subplot(1,2,1)
% Define the number of bins you want
numBins = 10; % Adjust this value as needed

% Create histograms for each group with reduced bins
histogram(mcineg_age, 'BinEdges', linspace(min([mcineg_age, mcipos_age])-10, max([mcineg_age, mcipos_age])+10, numBins+1), 'DisplayName', 'MCINeg', FaceColor=[0.9020    0.2941    0.2078]);
hold on;
histogram(mcipos_age, 'BinEdges', linspace(min([mcineg_age, mcipos_age])-10, max([mcineg_age, mcipos_age])+10, numBins), 'DisplayName', 'MCIPos', FaceColor=[0.3020    0.7333    0.8353]);
hold off;

xlabel('Age');
ylabel('Counts');
yticks([1,2,3,4,5])
legend;


[h,p,ci,stat] = ttest2(mcineg_age, mcipos_age);
% Create a string with t-test information
ttestInfo = sprintf('t(%d)= %.3f', stat.df, p);
% Add the t-test information to the figure
annotation('textbox', [0.15, 0.65, 0.1, 0.1], 'String', ttestInfo, 'BackgroundColor', 'w');

%exportgraphics(f,ResultFolder+"/AgeCompare.png",'Resolution',300);


%% compare the education duration between MCI+ and MCI-
subplot(1,2,2)
% Define the number of bins you want
numBins = 10; % Adjust this value as needed

% Create histograms for each group with reduced bins
histogram(mcineg_ed, 'BinEdges', linspace(min([mcineg_ed, mcineg_ed])-5, max([mcineg_ed, mcineg_ed])+5, numBins+1), 'DisplayName', 'MCINeg', FaceColor=[0.9020    0.2941    0.2078]);
hold on;
histogram(mcipos_ed, 'BinEdges', linspace(min([mcipos_ed, mcipos_ed])-5, max([mcipos_ed, mcipos_ed])+5, numBins), 'DisplayName', 'MCIPos', FaceColor=[0.3020    0.7333    0.8353]);
hold off;

xlabel('Years in education');
ylabel('Counts');
legend;

[h,p,ci,stat] = ttest2(mcineg_ed, mcipos_ed);
% Create a string with t-test information
ttestInfo = sprintf('t(%d)= %.3f', stat.df, p);
% Add the t-test information to the figure
annotation('textbox', [0.7, 0.65, 0.1, 0.1], 'String', ttestInfo, 'BackgroundColor', 'w');

exportgraphics(f,ResultFolder+"/AgeEducationCompare.png",'Resolution',300);

%% check the correlation between parameter values and the age or education
f = figure('Position', [100 100 600 600]);

colors = [0.9020    0.2941    0.2078; ...
          0.3020    0.7333    0.8353; ...
          0         0.6275    0.5294];
condlabel = ["no change", "no distal cue", "no optical flow"];
for cond = 1:3
    ho_g2_cond = HealthyControls.Results.estimatedParams{cond}(:,7);
    %remove NAN
    ho_g2_cond = ho_g2_cond(~isnan(ho_g2_cond));
    
    
    cc = corrcoef(ho_age, ho_g2_cond);
    corr_value = cc(1,2);
    
    scatter(ho_age, ho_g2_cond, 'o', 'filled', 'MarkerFaceColor',colors(cond,:))
    hold on;
    
    %Fit a linear model
    linear_fit = polyfit(ho_age, ho_g2_cond,1);
    y_fit = polyval(linear_fit, ho_age);
    
    %plot
    plot(ho_age, y_fit, 'r', 'LineWidth', 2, 'Color', colors(cond,:));
    
    xlabel('Age');
    ylabel('Parameter value');
%     % Display the correlation coefficient and linear fit equation on the plot
%     text(0.5, max(ho_g2_cond), ['Correlation: ' num2str(corr_value)], 'FontSize', 12, 'FontWeight', 'bold');
%     text(0.5, max(ho_g2_cond)-0.2, ['Linear Fit: y = ' num2str(linear_fit(1)) 'x + ' num2str(linear_fit(2))], 'FontSize', 12, 'FontWeight', 'bold');
%    

    %legend({'Data', 'Linear Fit'}, 'Location', 'best');
    
    
    % Calculate the number of data points
    n = length(ho_age);
    
    % Calculate the t-statistic
    t_statistic = corr_value * sqrt((n - 2) / (1 - corr_value^2));
    
    % Calculate the degrees of freedom for the t-distribution
    degrees_of_freedom = n - 2;
    
    % Calculate the p-value using the t-distribution
    p_value = 2 * (1 - tcdf(abs(t_statistic), degrees_of_freedom));
     
    fprintf('Correlation coefficient (r): %.2f, p-value: %.4f, t-statistic: %.4f\n', corr_value, p_value, t_statistic);
    %text(2.6, -0.1, condlabel(cond) + ' pvalue='+num2str(p_value), 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');

    hold on

end
    

