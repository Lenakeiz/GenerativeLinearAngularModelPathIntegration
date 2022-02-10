function PlotMentalTrajectory(Results, groupname,config)
%   Plot mental trajectory with estimated parameters across participants
%   ParamValues: G1, G2, G3, g2, g3, k3, sigma, nu
%   X: the Actual positions
%   DX:
%   THETADX: 

%extract data
ParamValues = Results.estimatedParams;
X = Results.X;
DX = Results.DX;
THETADX = Results.THETADX;
flagOoB = Results.flagOoB;

%load configurations necessary for the script
trialfolder = config.TrialFolder;

disp("%%%%%%%%%%%%%%% PLOT MENTAL TRAJECTORY OF GROUP: " + groupname+ "%%%%%%%%%%%%%%%");

%create storing folder for trajectory if not exist
groupfoldername = trialfolder+groupname+'/Mental/';
if ~exist(groupfoldername, 'dir')
   mkdir(groupfoldername);
end

subjectSize = size(X,2);

for subj = 1:subjectSize
    
    disp("subject id is "+num2str(subj))
    subjX = X{subj};
    subjDX = DX{subj};
    subjTHETADX = THETADX{subj};
    paramX = ParamValues(subj,:);
    subjflagOoB = flagOoB{subj};

    [G1,G2,G3,g2,g3,k3,sigma, nu] = deal(paramX(1),paramX(2),paramX(3),paramX(4),paramX(5),paramX(6),paramX(7),paramX(8)); 

    %for each subjects, lets randomly sample 8 trials for visualizing
    sampleSize = size(subjX,2);
    %idxs = randsample(sampleSize, 8);
    %for each subjects, lets visualize the first 8 trials
    %idxs = 1:8;
    
    %generate a empty figure for storing data
    f = figure('visible','off','Position', [100 100 1000 800]);
    
    plot_row = ceil(sqrt(sampleSize));
    plot_col = ceil(sqrt(sampleSize));
    
    %tiledlayout(plot_row, plot_col, 'TileSpacing', 'none', 'Padding', 'none');

    for tr_id=1:sampleSize
        
        subplot(plot_row,plot_col,tr_id);
        %nexttile 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %first draw the physical trajectory%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x_ = subjX{tr_id}(:,1);
        y_ = subjX{tr_id}(:,2);
        plot(x_', y_'-0.1, '.-', 'Markersize', 10, LineWidth=1, Color='k');
        hold on

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %then draw the mental trajectory%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %leg length
        l1 = subjDX{tr_id}(1); 
        l2 = subjDX{tr_id}(2);
        %angle value
        theta2 = subjTHETADX{tr_id}(2);
        theta2_prime = g2*theta2;

        %mental leg 1 start point:
        start = [0,0];
        %mental leg 1 end point:
        men_p1 = [G1*l1, 0];
        %mental leg 2 end point:
        men_p2 = [G1*l1+G2*l2*cos(theta2_prime), G2*l2*sin(theta2_prime)];

        %for mental leg3 end point:

        %first calculate length of accurate mental vector 3
        h = norm(men_p2);
        %second calculate turn angle of accurate mental vector 3
        vec1 = men_p2-men_p1; vec2 = [0,0]-men_p2;
        alpha = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
        %transfer from degree to radians
        alpha = deg2rad(alpha);
        %wrap to (0,2pi)
        alpha = mod(alpha, 2*pi);

        %third considering  execution turn error
        theta3_prime = g3*alpha+k3*(1-g3);
        %also wrap to (0,2pi)
        theta3_prime = mod(theta3_prime, 2*pi);       

        %mental vector 3 with execution error
        exe_men_vec3 = [G3*h*cos(theta2_prime+theta3_prime), G3*h*sin(theta2_prime+theta3_prime)];
        men_p3 = men_p2 + exe_men_vec3;

        %draw mental trajectory
        mentalxy = [start;men_p1;men_p2;men_p3];
        m_x = mentalxy(:,1); m_y = mentalxy(:,2);
		
        plot(m_x', m_y', '.-', 'Markersize', 10,  LineWidth=1, Color='r');

        hold on

        %draw mental trajectory start from the physical end point 3
        mental_vec3_from_phyiscal = [G3*h*cos(theta2+theta3_prime), G3*h*sin(theta2+theta3_prime)];
        start = [x_(3),y_(3)];
        men_p3_from_physical = start+mental_vec3_from_phyiscal;
        mentalxy_from_physical = [start;men_p3_from_physical];
        m_x_from_physical = mentalxy_from_physical(:,1); 
        m_y_from_physical = mentalxy_from_physical(:,2);
        plot(m_x_from_physical', m_y_from_physical', '.-', 'Markersize', 10,  LineWidth=1, Color='g');

        hold on

%         %a multivarite normal distirbution around zero point with estimated sigma
%         mu = [0,0]; Sigma = [sigma, 0; 0, sigma]; %this covariance matrix will be estimated in more detailed way in the future 
%         %lets say, generate 1000 samples
%         sampling_dots = mvnrnd(mu, Sigma, 1000);
% 
%         %scatter plot of the sampling
%         scatter(sampling_dots(:,1), sampling_dots(:,2), 2, 'b', 'filled', 'MarkerFaceColor','b','MarkerEdgeColor','b', 'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',0.05);
        
        if subjflagOoB(tr_id)==0 %non-oob trial
            title("tr="+num2str(tr_id)+' \alpha='+num2str(round(alpha,2))+' \theta_3^\prime='+num2str(round(theta3_prime,2)));
        else %oob trial
            title("\color{red}tr="+num2str(tr_id)+' \alpha='+num2str(round(alpha,2))+' \theta_3^\prime='+num2str(round(theta3_prime,2)));
        end
        %xlabel('X position'); ylabel('Y position');

        axis equal
        xlim([-6,6]);ylim([-6,6]);

        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        box off

        %set axis with equal scale
        %set(gca,'DataAspectRatio',[1 1 1]);
    end

    sgtitle("Group name: "+ groupname + ...
            " Subject ID: "+num2str(subj)+ ...
            ' G_1='+num2str(round(G1,2))+ ...
            ' G_2='+num2str(round(G2,2))+ ...
            " g_2="+num2str(round(g2,2))+ ...
            " g_3="+num2str(round(g3,2))+ ...
            " k_3="+num2str(round(k3,2)), 'FontSize', 20);

    exportgraphics(f,groupfoldername+string(subj)+".png",'Resolution',300);

    close(f);

end

end
