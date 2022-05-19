function VisualizeRealtimeTrackingData(Group, pId, trialId, playbackSpeed,varargin)

% Reconstructing a visualization of participants walk and head drection after reaching the third
% cone
config.cutFromConeThree = false;

if exist('varargin','var')
    for i = 1:2:nargin-4
        if(strcmpi(varargin{i},'cutconethree'))
            config.cutFromConeThree = varargin{i+1};
        end
    end
end

ColorPattern;

Cone_pos      = Group.FlagPos{1,pId}{trialId,1};
Extracted_pos = Group.Path{1,pId}{trialId,1};
Extracted_pos = array2table(Extracted_pos,"VariableNames",{'Time' 'Pos_X' 'Pos_Y' 'Pos_Z' 'Forward_X' 'Forward_Y' 'Forward_Z'});

close all; clc;

% Extracting after reaching cone 3
if(config.cutFromConeThree)
    Extracted_pos = Extracted_pos(Extracted_pos.Time > Group.FlagTrigTimes{1,pId}{trialId,3},:);
end

% Plotting according to the speed
CreateCustomFigure;

hold all

xlim([-2.5 2.5])
ylim([-2.5 2.5])

xlabel('x (m)')
ylabel('y (m)')

axis square

tracking_size = height(Extracted_pos);

% Plot the position of the third cone
if(config.cutFromConeThree == false)
    plot(Cone_pos(1,1),Cone_pos(1,3),'Marker','d','LineStyle','none','MarkerFaceColor',config.color_scheme_npg(2,:),'MarkerEdgeColor','black','MarkerSize',12); 
    plot(Cone_pos(2,1),Cone_pos(2,3),'Marker','d','LineStyle','none','MarkerFaceColor',config.color_scheme_npg(3,:),'MarkerEdgeColor','black','MarkerSize',12); 
end
plot(Cone_pos(3,1),Cone_pos(3,3),'Marker','d','LineStyle','none','MarkerFaceColor',config.color_scheme_npg(1,:),'MarkerEdgeColor','black','MarkerSize',12); 

for k = 2 : tracking_size

    tempPos = Extracted_pos(k-1:k,[2 4]);
    
    % Current path
    plot(tempPos.Pos_X,tempPos.Pos_Z,'Marker','none','LineStyle','-','LineWidth',0.7,'Color',config.color_scheme_npg(4,:));
    
    % Adding an arrow to visualize the orientation of the participant at
    % the second point    
    % Extracting direction on the xz plane (the projection on the floor)
    tempDir = [Extracted_pos.Forward_X(k) Extracted_pos.Forward_Z(k)];
    tempDir = tempDir/norm(tempDir);
    % Making it smaller for visualization
    tempDir = tempDir.*0.2;

    qv = quiver(Extracted_pos.Pos_X(k),Extracted_pos.Pos_Z(k),tempDir(1), tempDir(2), 'off');
    qv.Color = config.color_scheme_npg(5,:);
    qv.LineWidth = 1.3;
    
    pause(playbackSpeed);

    delete(qv);

end

hold off

end