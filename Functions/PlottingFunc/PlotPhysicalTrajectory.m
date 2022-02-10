function PlotPhysicalTrajectory(Results, groupname,config)
%   Plot mental trajectory with estimated parameters across participants
%   X: the Actual positions

%extract data
X = Results.X;
flagOoB = Results.flagOoB;

%load configurations necessary for the script
trialfolder = config.TrialFolder;

disp("%%%%%%%%%%%%%%% PLOT PHYSICAL TRAJECTORY OF GROUP: " + groupname+ "%%%%%%%%%%%%%%%");

%create storing folder for trajectory if not exist
groupfoldername = trialfolder+groupname+'/Physical/';
if ~exist(groupfoldername, 'dir')
   mkdir(groupfoldername);
end

subjectSize = size(X,2);

for subj = 1:subjectSize
    
    disp("subject id is "+num2str(subj))
    subjX = X{subj};
    subjflagOoB = flagOoB{subj};

    sampleSize = size(subjX,2);
    
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
        axis equal
        xlim([-6,6]);ylim([-6,6]);

        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        box off

        if subjflagOoB(tr_id)==0 %non-oob trial
            title("tr="+num2str(tr_id));
        else %oob trial
            title("\color{red}tr="+num2str(tr_id));
        end        
    end

    sgtitle("Group name: "+ groupname +"    Subject ID: "+num2str(subj), 'FontSize', 20);

    exportgraphics(f,groupfoldername+string(subj)+".png",'Resolution',300);

    close(f);

end

end
