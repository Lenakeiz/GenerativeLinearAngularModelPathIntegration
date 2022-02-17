function PlotEndpointDistribution(Results, groupname, config)
%   Error bar plot of estimated parameters across participants
%   ParamValues: beta1, beta2, beta3, thetaB, thetaD, sigma
%   X: the Actual positions
%   DX:
%   THETADX: 


%extract data
ParamValues = Results.estimatedParams;
X = Results.X;
DX = Results.DX;
THETADX = Results.THETADX;

%load configurations necessary for the script
trialfolder = config.TrialFolder;

%sample 8 subjects in the group
subjectSize = size(X,2);
%idxs = randsample(subjectSize, 8);

%visualize the first 8 subjects
%idxs = 1:8;

%generate a empty figure for storing data
f = figure('visible','off','Position', [100 100 2000 1000]);

plot_row = ceil(sqrt(subjectSize));
plot_col = 2*ceil(sqrt(subjectSize));

for subj_id = 1:subjectSize

    disp("subject id is "+num2str(subj_id))
    subjX = X{subj_id};
    subjDX = DX{subj_id};
    subjTHETADX = THETADX{subj_id};
    paramX = ParamValues(subj_id,:);

    [G1,G2,G3,g2,g3,k3,sigma, nu] = deal(paramX(1),paramX(2),paramX(3),paramX(4),paramX(5),paramX(6),paramX(7),paramX(8));

    %for each subjects, calculate the number of trial
    sampleSize = size(subjX,2);
    
    Phys_EndPoints = zeros(sampleSize,2);
    Ment_EndPoints = zeros(sampleSize,2);

    for tr_id=1:sampleSize

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %extract the physical end point%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        physical_endpoint = subjX{tr_id}(end,:);
        Phys_EndPoints(tr_id,:) = physical_endpoint;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %calculate the mental end point%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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

        %execution mental vector 3
        exe_men_vec3 = [G3*h*cos(theta2_prime+theta3_prime), G3*h*sin(theta2_prime+theta3_prime)];
        men_p3 = men_p2 + exe_men_vec3;

        %store the end point
        Ment_EndPoints(tr_id,:) = men_p3;
    end

    subplot(plot_row,plot_col,2*subj_id-1);  
    %plot all physical end point 
    scatter(Phys_EndPoints(:,1), Phys_EndPoints(:,2), 5, 'filled', 'MarkerFaceColor','b','MarkerEdgeColor','b', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0.5);
    hold on;
    %plot all mental end point 
    scatter(Ment_EndPoints(:,1), Ment_EndPoints(:,2), 5, 'filled', 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',0.5);        
    hold on
    
%     %plot a circle with one sigma
%     r = 0:0.01:2*pi;
%     x = sigma*cos(r); y = sigma*sin(r);
%     scatter(x, y, 1, 'MarkerFaceColor','k','MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)

    axis equal

    title("ID: "+num2str(subj_id), 'FontSize', 8);
    xtick = get(gca,'XTickLabel');
    set(gca,'XTickLabel',xtick, 'fontsize',6);
    ytick = get(gca,'YTickLabel');
    set(gca,'YTickLabel',ytick, 'fontsize',6);       
    hold off;
    
    %plot the corresponding estimated parameters
    subplot(plot_row,plot_col,2*subj_id);            
    %draw mental trajectory with "mul" with regression to k
    HX = categorical({'G_1','G_2','G_3','g_2','g_3','k_3'});
    HX = reordercats(HX,{'G_1','G_2','G_3','g_2','g_3','k_3'});
    HY = [round(G1,2), round(G2,2), round(G3,2), round(g2,2), round(g3,2), round(k3,2)];

    b = bar(HX,HY);
    xtips = b.XEndPoints;
    ytips = b.YEndPoints;
    labels = string(b.YData);
    text(xtips,ytips,labels,'HorizontalAlignment','center','VerticalAlignment','bottom', 'FontSize',6);
    hold on;
    yline(1);
    hold off;
    set(gca,'box','off');
    xtick = get(gca,'XTickLabel');
    set(gca,'XTickLabel',xtick, 'fontsize',6);
    ytick = get(gca,'YTickLabel');
    set(gca,'YTickLabel',ytick, 'fontsize',6);         
    %set axis with equal scale
    %set(gca,'DataAspectRatio',[1 1 1]);   
end

sgtitle("Group name: "+ groupname);

%extract figure with high resolution

figname = trialfolder+"EndPointDistribution_"+groupname + ".png";
exportgraphics(f,figname,'Resolution',300);

close(f);

end
