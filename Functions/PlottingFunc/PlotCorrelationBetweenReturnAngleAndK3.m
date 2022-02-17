function AllTurnAngleDiff=PlotCorrelationBetweenReturnAngleAndK3(AllParamValues, AllDX, AllTHETADX, groupname, config)
%%plot last turn distribution
% ParamValues: estimated parameter values: numSubjects x numParams
% DX is a cell structure containing the segment of each trial
% THETAX is the turning angle (wrong at the moment)
% groupname: self-explaining
% config: model configuration

resultfolder = config.ResultFolder;

%create storing folder for trajectory if not exist
savefoldername = resultfolder+"CorrelationBetweenReturnAngleAndK3/";
if ~exist(savefoldername, 'dir')
   mkdir(savefoldername);
end

%generate a empty figure for storing data
f = figure('visible','off','Position', [100 100 1000 600]);
trialname = ["No change", "No distal cue", "No optical flow"];

numConds = length(AllParamValues);

AllTurnAngleDiff = cell(0);

for TRIAL_FILTER=1:numConds
    
    ParamValues = AllParamValues{TRIAL_FILTER};
    DX = AllDX{TRIAL_FILTER};
    THETADX = AllTHETADX{TRIAL_FILTER};
    
    MeanTurnAngle = []; k3_all = []; g3_all = [];
    
    subjectSize = size(THETADX,2);

    subplot(2,3,TRIAL_FILTER)
    
    for subj=1:subjectSize
        paramX = ParamValues(subj,:);
        subjTHETADX = THETADX{subj};
        subjDX = DX{subj};
    
        %extract parameters
        [gamma,G3,g2,g3,k3,sigma,nu] = deal(paramX(1),paramX(2),paramX(3),paramX(4),paramX(5),paramX(6),paramX(7));
        %k3 = k3/(1+eps-g3);
        k3_all = [k3_all, k3];
        g3_all = [g3_all, g3];
        
        %get the turning angle from all trials from one subject
        sampleSize = size(subjTHETADX,2);
        subjtheta3 = zeros(1,sampleSize);
        repeat_k3 = zeros(1,sampleSize);
        for tr_id=1:sampleSize
            theta3 = subjTHETADX{tr_id}(3);
            repeat_k3(tr_id) = k3;
            subjtheta3(tr_id) = theta3;
        end
        MeanTurnAngle = [MeanTurnAngle,mean(subjtheta3)];
        %plot the scatter of return angle of a subject
        scatter(subjtheta3+0.2, repeat_k3, 10, 'dr', 'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
        hold on
    end
    
    hold on
    %do linear regression 
    md1 = fitlm(MeanTurnAngle,k3_all, 'linear'); %returns a linear regression model of the responses K3, fit to MeanTurnAngle.
    h = plot(md1, 'LineWidth', 3);
    set(h(1), 'Color', 'r', 'Marker', '.', 'MarkerSize', 20); %set data point
    set(h(2), 'Color', 'k', 'LineWidth', 3); %set the regression line
    set(h(3), 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.'); %set the upper confidence line
    set(h(4), 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.'); %set the lower confidence line
    %hold on
    %scatter(MeanTurnAngle, k3_all, 20, 'dg', 'filled');
    %pearson correlation
    [R,P]=corrcoef(MeanTurnAngle,k3_all);
    %text(0.1,5.5, "corrcoef="+num2str(round(R(1,2),2)), 'FontSize', 10, 'Color', 'r');
    %text(0.1,5.0, "pvalue="+num2str(round(P(1,2),2)), 'FontSize', 10, 'Color', 'r');

    axis equal
    legend('off');
    %legend('Location','northwest', 'FontSize',5);
    %xlim([0,2*pi]);ylim([0,2*pi]);
    ylabel('Preferred angle $k_3$', 'Interpreter','latex', 'FontSize', 10);
    xlabel('Mean return angle $\bar{\theta}_3$', 'Interpreter','latex', 'FontSize', 10);
    title(trialname(TRIAL_FILTER), 'FontSize', 15);

    subplot(2,3,3+TRIAL_FILTER)
    axis equal
    box on
    scatter(k3_all, g3_all, 20, 'ob', 'filled');
    %xlim([0,2*pi]);ylim([0,1]);
    %xticks([0,pi,2*pi]);
    %xticklabels(["0" "\pi" "2\pi"]);
    xlabel('Preferred angle $k_3$', 'Interpreter','latex', 'FontSize', 10);
    ylabel('Estimated $g_3$', 'Interpreter','latex', 'FontSize', 10);   

    AllTurnAngleDiff{TRIAL_FILTER}=k3_all-MeanTurnAngle;
end

sgtitle("GroupName: "+groupname, 'FontSize', 20);
    
exportgraphics(f,savefoldername+ groupname+".png",'Resolution',300);

end