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

ax = gca;
ax.FontSize = 20;
xlim([-2.5 2.5])
ylim([-2.5 2.5])

xlabel('x (m)')
ylabel('y (m)')

if(isfield(Group,'Reconstructed'))
    calculatedAngle = Group.Reconstructed{1,pId}.InboundRotation(trialId);
    title(['Participant: ', num2str(pId), ' Trial: ', num2str(trialId), ' Body Rotation: ', num2str(calculatedAngle,'%.0f') ],FontSize=25);
else
    title(['Participant: ', num2str(pId), ' Trial ', num2str(trialId)],FontSize=25);
end

axis square

tracking_size = height(Extracted_pos);

% Plot the position of the third cone
pCone1 = plot(Cone_pos(1,1),Cone_pos(1,3),'Marker','d','LineStyle','none','MarkerFaceColor',config.color_scheme_npg(2,:),'MarkerEdgeColor','black','MarkerSize',18); 
pCone2 = plot(Cone_pos(2,1),Cone_pos(2,3),'Marker','d','LineStyle','none','MarkerFaceColor',config.color_scheme_npg(3,:),'MarkerEdgeColor','black','MarkerSize',18); 
pCone3 = plot(Cone_pos(3,1),Cone_pos(3,3),'Marker','d','LineStyle','none','MarkerFaceColor',config.color_scheme_npg(1,:),'MarkerEdgeColor','black','MarkerSize',18); 

leg = legend([pCone1 pCone2 pCone3], 'Cone 1', 'Cone 2', 'Cone 3', 'AutoUpdate','off');
leg.Location = 'northeastoutside';
leg.FontSize = 15;

for k = 2 : tracking_size

    tempPos = Extracted_pos(k-1:k,[2 4]);
    
    % Current path
    plot(tempPos.Pos_X,tempPos.Pos_Z,'Marker','none','LineStyle','-','LineWidth',2.5,'Color',config.color_scheme_npg(4,:));
    
    % Adding an arrow to visualize the orientation of the participant at
    % the second point    
    % Extracting direction on the xz plane (the projection on the floor)
    tempDir = [Extracted_pos.Forward_X(k) Extracted_pos.Forward_Z(k)];
    tempDir = tempDir/norm(tempDir);
    % Making it smaller for visualization
    tempDir = tempDir.*0.5;

    qv = quiver(Extracted_pos.Pos_X(k),Extracted_pos.Pos_Z(k),tempDir(1), tempDir(2), 'off');
    qv.Color = config.color_scheme_npg(5,:);
    qv.LineWidth = 5;
    qv.MaxHeadSize = 2;
    
    drawnow

    pause(playbackSpeed);

    if(k<tracking_size)
        delete(qv);
    end

end

hold off

end