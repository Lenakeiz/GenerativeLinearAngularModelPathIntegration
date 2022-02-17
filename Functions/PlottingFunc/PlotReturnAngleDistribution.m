function PlotReturnAngleDistribution(AllParamValues, AllDX, AllTHETADX, AllflagOoB, groupname, config)
%%plot return angle distribution of theta3, alpha, theta3prime, and sampling results

% ParamValues: estimated parameter values: numSubjects x numParams
% DX is a cell structure containing the segment of each trial
% THETAX is the turning angle (wrong at the moment)
% flagOoB is the flag of Out-of-Boundary trials, with 0 non-OoB trials 1 those OoB trials
% groupname: self-explaining
% config: model configuration

resultfolder = config.ResultFolder;

%create storing folder for trajectory if not exist
savefoldername = resultfolder+"ReturnAngleDistribution/";
if ~exist(savefoldername, 'dir')
   mkdir(savefoldername);
end

trialname = ["No change", "No distal cue", "No optical flow"];

numConds = length(AllParamValues);

for TRIAL_FILTER=1:numConds
    
    ParamValues = AllParamValues{TRIAL_FILTER};
    DX = AllDX{TRIAL_FILTER};
    THETADX = AllTHETADX{TRIAL_FILTER};
    flagOoB = AllflagOoB{TRIAL_FILTER};
    
    samp_num = 10;

    subjectSize = size(THETADX,2);

    all_angles = []; all_alphas = []; all_theta_primes = []; all_flagOoBs = [];
    all_samples = [];
    all_g3 = []; all_b = [];
    total_angle = cell(0); total_alpha = cell(0); 

    %generate a empty figure for storing data for each subject
    f = figure('visible','off','Position', [100 100 2000 1000]);
    
    plot_row = ceil(sqrt(subjectSize));
    plot_col = ceil(sqrt(subjectSize));

    for subj=1:subjectSize
        paramX = ParamValues(subj,:);
        subjTHETADX = THETADX{subj};
        subjDX = DX{subj};
        subjflagOoB = flagOoB{subj};
    
        %extract parameters
        [gamma,G3,g2,g3,b,sigma,nu] = deal(paramX(1),paramX(2),paramX(3),paramX(4),paramX(5),paramX(6),paramX(7));
        all_g3 = [all_g3, g3]; all_b = [all_b, b];
        G1 = gamma^2*G3;
        G2 = gamma*G3;  
    
        sampleSize = size(subjTHETADX,2);
        Angles = zeros(1,sampleSize);
        Alphas = zeros(1,sampleSize); 
        Theta_primes = zeros(1,sampleSize);
        Samples = zeros(1, sampleSize*samp_num);

        %get theta3_bar
        theta3_bar_all = zeros(1,sampleSize);
        for tr = 1:sampleSize
            theta3 = subjTHETADX{tr}(3);
            theta3_bar_all(tr) = theta3;
        end
        theta3_bar = mean(theta3_bar_all);
    
        for tr_id=1:sampleSize
    
            turn3_angle = subjTHETADX{tr_id}(3);
            Angles(tr_id) = turn3_angle;
    
            %accurate mental angle, i.e., alpha
            %leg length
            l1 = subjDX{tr_id}(1); 
            l2 = subjDX{tr_id}(2);
            %angle value
            theta2 = subjTHETADX{tr_id}(2);
            theta2_prime = g2*theta2;
        
            %mental leg 1 end point:
            men_p1 = [G1*l1, 0];

            %mental leg 2 end point:
            men_p2 = [G1*l1+G2*l2*cos(theta2_prime), G2*l2*sin(theta2_prime)];
        
            %for mental leg3 end point:
            %calculate turn angle of accurate mental vector 3
            vec1 = men_p2-men_p1; vec2 = [0,0]-men_p2;
            alpha = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
            %transfer from degree to radians
            alpha = deg2rad(alpha);
            %wrap to (0,2pi)
            alpha = mod(alpha, 2*pi);  
    
            Alphas(tr_id)=alpha;
    
            %considering execution turn error
            theta3_prime = g3*alpha+b;
            %also wrap to (0,2pi)
            %theta3_prime = mod(theta3_prime, 2*pi);     
    
            Theta_primes(tr_id) = theta3_prime;
    
            %considering sampling, for each theta3_prime, sampling 'samp_num' times
            %according to the von mises distribution nu
            sampled_theta = circ_vmrnd(theta3_prime, nu, samp_num); %column vector
            sampled_theta = sampled_theta';
            %mod to (0,2*pi)
            sampled_theta = mod(sampled_theta, 2*pi);
            Samples((tr_id-1)*samp_num+1:tr_id*samp_num)=sampled_theta;
        end

        all_angles = [all_angles,Angles];
        all_alphas = [all_alphas, Alphas];
        all_theta_primes = [all_theta_primes, Theta_primes];
        all_flagOoBs = [all_flagOoBs, subjflagOoB]; all_flagOoBs = logical(all_flagOoBs); %turn to logcial value for indexing
        all_samples = [all_samples, Samples];
        total_angle{subj} = Angles;
        total_alpha{subj} = Alphas;

        %each subject will have a sub figure
        subplot(plot_row,plot_col,subj);
    
        %for each participants, plot theta_3_prime against alpha 
        scatter(Alphas, Theta_primes, 2, 'b', 'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0.5);
        hold on
        scatter(Alphas, Angles, 2, 'r', 'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0.5);
        hold on
        %add refrence line of y=1*x+0
        plot([0,2*pi],[0,2*pi], 'Color','k','LineStyle','--', 'LineWidth',0.5);
        xlabel('Encoded \alpha'); ylabel('Executed \theta_3^\prime or physical \theta_3');
        axis equal
        xticks([ 0, pi, 2*pi]);yticks([ 0, pi, 2*pi]);
        xticklabels({'0','\pi','2\pi'});yticklabels({'0','\pi','2\pi'});
        xlim([0,2*pi]);ylim([0,2*pi]);
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        box off
        if g3<0.05
            title("\color{red}id=" + num2str(subj) + ...
                  " G_1="+num2str(round(G1,2))+...
                  " G_2="+num2str(round(G2,2))+...
                  " g_3="+num2str(round(g3,2))+ ...
                  " \theta="+num2str(round(theta3_bar,2))+...
                  " b="+num2str(round(b,2)),'FontSize',6);
        else
            title("id=" + num2str(subj) + ...
                  " G_1="+num2str(round(G1,2))+...
                  " G_2="+num2str(round(G2,2))+...
                  " g_3="+num2str(round(g3,2))+ ...
                  " \theta="+num2str(round(theta3_bar,2))+...
                  " b="+num2str(round(b,2)),'FontSize',6);
        end
    end
                                                        

    sgtitle(groupname +" "+trialname(TRIAL_FILTER) +' Mean g_3='+num2str(round(mean(all_g3),2))+' Mean c='+num2str(round(mean(all_b),2)));
    exportgraphics(f,savefoldername+groupname+"_"+trialname(TRIAL_FILTER)+"_AngleRegression.png",'Resolution',300);
    close(f);

    %% plot all subjects' thets3 against alpha in different colors  
    f = figure('visible','off','Position', [100 100 500 500]);
    cm = jet(subjectSize); %pick numSubjs colors from 'jet'
    for subj=1:subjectSize
        scatter(total_alpha{subj}, total_angle{subj}, 20, cm(subj,:), 'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0.5)
        hold on
    end
    %add refrence line of y=1*x+0
    plot([0,2*pi],[0,2*pi], 'Color','k','LineStyle','--', 'LineWidth',3);  
    hold on

    %do linear regression 
    md1 = fitlm(all_alphas,all_angles, 'linear'); %returns a linear regression model of the responses all_angles, fit to all_alphas.
    h = plot(md1, 'LineWidth', 3);
    set(h(1), 'Color', 'k', 'Marker', '.', 'MarkerSize', 0.1); %set data point
    set(h(2), 'Color', 'r', 'LineWidth', 2); %set the regression line
    set(h(3), 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.'); %set the upper confidence line
    set(h(4), 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.'); %set the lower confidence line
    
    P = polyfit(all_alphas,all_angles,1);
    slope = P(1);
    %text(0.1,5.5, "slope="+num2str(round(slope,4)), 'FontSize', 10, 'Color', 'r');
    
    legend('off');
    xlabel("Encoded \alpha"); ylabel("physical \theta_3");
    axis equal
    xticks([ 0, pi, 2*pi]);yticks([ 0, pi, 2*pi]);
    xticklabels({'0','\pi','2\pi'});yticklabels({'0','\pi','2\pi'});
    xlim([0,2*pi]);ylim([0,2*pi]);
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off   

    title(groupname +" "+trialname(TRIAL_FILTER)+" Slope="+num2str(round(slope,4)));
    exportgraphics(f, savefoldername+groupname+"_"+trialname(TRIAL_FILTER)+"__AngleRegressionAll.png",'Resolution',300)
    close(f);

    %% plot return angle distribution
    f = figure('visible','off','Position', [100 100 1000 1000]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot histogram of all subjects' third turn in subplots
    sgtitle("Group name: "+ groupname + " "+trialname(TRIAL_FILTER));
    
    subplot(4,1,1)
    %non OoB trials
    histogram(all_angles(~all_flagOoBs), 30, 'FaceColor','m', 'facealpha',0.5);
    hold on
    %OoB trials
    histogram(all_angles(all_flagOoBs), 30, 'FaceColor','k', 'facealpha',0.5);
    ylabel('frequency');xlabel('physical angle', 'FontSize', 20);
    xticks([ 0, pi/2, pi, 3/2*pi, 2*pi]);
    xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
    xlim([0, 2*pi]);
    
    subplot(4,1,2)
    histogram(all_alphas(~all_flagOoBs), 30, 'FaceColor','b', 'facealpha',0.5);
    hold on
    histogram(all_alphas(all_flagOoBs), 30, 'FaceColor','k', 'facealpha',0.5);
    ylabel('frequency');xlabel('accurate mental angle', 'FontSize', 20);
    xticks([ 0, pi/2, pi, 3/2*pi, 2*pi]);
    xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
    xlim([0, 2*pi])
    
    subplot(4,1,3)
    histogram(all_theta_primes(~all_flagOoBs), 30, 'FaceColor','g', 'facealpha',0.5);
    hold on
    histogram(all_theta_primes(all_flagOoBs), 30, 'FaceColor','k', 'facealpha',0.5);
    ylabel('frequency');
    xlabel('mental angle with execution error', 'FontSize', 20); 
    xticks([ 0, pi/2, pi, 3/2*pi, 2*pi]);
    xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
    xlim([0, 2*pi])
    
    subplot(4,1,4)
    histogram(all_samples, 30, 'FaceColor','c', 'facealpha',0.5);
    ylabel('frequency');
    xlabel('sampling from a Von Mises distribution with mean \theta_3^\prime and variance \nu', 'FontSize', 20); 
    xticks([ 0, pi/2, pi, 3/2*pi, 2*pi]);
    xticklabels({'0','\pi/2','\pi', '3/2\pi','2\pi'});
    xlim([0, 2*pi])

    exportgraphics(f,savefoldername+groupname+"_"+trialname(TRIAL_FILTER)+".png",'Resolution',300);
    close(f);
end



end