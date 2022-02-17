function PlotReturnAngleDistributionInVerticalScatter(AllParamValues, AllDX, AllTHETADX, groupname, config)
%%plot last turn distribution
% ParamValues: estimated parameter values: numSubjects x numParams
% DX is a cell structure containing the segment of each trial
% THETAX is the turning angle (wrong at the moment)
% groupname: self-explaining
% config: model configuration

resultfolder = config.ResultFolder;

%create storing folder for trajectory if not exist
savefoldername = resultfolder+"ReturnAngleDistributionInVerticalScatter/";
if ~exist(savefoldername, 'dir')
   mkdir(savefoldername);
end

%generate a empty figure for storing data
f = figure('visible','off','Position', [100 100 1000 400]);
trialname = ["No change", "No distal cue", "No optical flow"];

numConds = length(AllParamValues);

for TRIAL_FILTER=1:numConds

    subplot(1,3,TRIAL_FILTER)
    
    ParamValues = AllParamValues{TRIAL_FILTER};
    DX = AllDX{TRIAL_FILTER};
    THETADX = AllTHETADX{TRIAL_FILTER};
    
    all_theta3 = []; all_alpha = []; all_theta3_primes = []; all_c = [];
    
    subjectSize = size(THETADX,2);
    cm = jet(subjectSize); %pick numSubjs colors from 'jet'
    
    for subj=1:subjectSize
        paramX = ParamValues(subj,:);
        subjTHETADX = THETADX{subj};
        subjDX = DX{subj};
    
        %extract parameters
        [gamma,G3,g2,g3,c,sigma,nu] = deal(paramX(1),paramX(2),paramX(3),paramX(4),paramX(5),paramX(6),paramX(7));
        all_c = [all_c, c];
    
        sampleSize = size(subjTHETADX,2);
        Theta3 = zeros(1,sampleSize);
        Alpha = zeros(1, sampleSize);
        Theta3_primes = zeros(1,sampleSize);

        %get theta3_bar
        theta3_bar_all = zeros(1,sampleSize);
        for tr = 1:sampleSize
            theta3 = subjTHETADX{tr}(3);
            theta3_bar_all(tr) = theta3;
        end
        theta3_bar = mean(theta3_bar_all);

        for tr_id=1:sampleSize
    
            theta3 = subjTHETADX{tr_id}(3);
            Theta3(tr_id) = theta3;
    
            %accurate mental angle, i.e., alpha
            %leg length
            l1 = subjDX{tr_id}(1); 
            l2 = subjDX{tr_id}(2);
            %angle value
            theta2 = subjTHETADX{tr_id}(2);
            theta2_prime = g2*theta2;
        
            %mental leg 1 end point:
            G1 = gamma^2*G3;
            men_p1 = [G1*l1, 0];

            %mental leg 2 end point:
            G2 = gamma*G3; 
            men_p2 = [G1*l1+G2*l2*cos(theta2_prime), G2*l2*sin(theta2_prime)];
        
            %for mental leg3 end point:
            %calculate turn angle of accurate mental vector 3
            vec1 = men_p2-men_p1; vec2 = [0,0]-men_p2;
            alpha = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
            %transfer from degree to radians
            alpha = deg2rad(alpha);
            %wrap to (0,2pi)
            alpha = mod(alpha, 2*pi);
            Alpha(tr_id)=alpha;
    
            %considering execution turn error
            theta3_prime = g3*alpha+(1-g3)*theta3_bar+c;
            %also wrap to (0,2pi)
            theta3_prime = mod(theta3_prime, 2*pi);     
    
            Theta3_primes(tr_id) = theta3_prime;

            pt = plot([1, 2, 3], [theta3, alpha, theta3_prime], '-o', "Color", cm(subj,:), 'linewidth',2, 'MarkerSize', 0.5);
            %transparancy
            pt.Color = [pt.Color 0.05];
            hold on
        end
        all_theta3 = [all_theta3,Theta3];
        all_alpha = [all_alpha, Alpha];
        all_theta3_primes = [all_theta3_primes, Theta3_primes];
    end

    theta3_mean = mean(all_theta3);
    theta3_prime_mean = mean(all_theta3_primes);
    alpha_mean = mean(all_alpha);

    plot([1 2 3], [theta3_mean alpha_mean theta3_prime_mean], '-dk', 'linewidth',5, 'MarkerSize', 5)

    ylabel('Anguar value', 'FontSize', 15);
    xlim([0.5,4.5]);
    ylim([0,2*pi]);
    xticks([1 2 3 4]);
    xticklabels({'\theta_3','\alpha', 'k_3', '\theta_3^\prime', });
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',10);
    title(trialname(TRIAL_FILTER), 'FontSize', 15)

end

sgtitle("GroupName: "+groupname, 'FontSize', 20);
    
exportgraphics(f,savefoldername+ groupname+".png",'Resolution',300);

end