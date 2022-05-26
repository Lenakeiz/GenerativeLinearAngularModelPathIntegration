function VisualizeRealtimeTrackingData(Group, pId, trialId, playbackSpeed,varargin)

anglebetween = @(v,w) atan2d(w(:,2).*v(:,1) - v(:,2).*w(:,1), v(:,1).*w(:,1) + v(:,2).*w(:,2));

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
Trig_pos      = Group.TrigPos{1,pId}{trialId,1};
Extracted_pos = Group.Path{1,pId}{trialId,1};
Extracted_pos = array2table(Extracted_pos,"VariableNames",{'Time' 'Pos_X' 'Pos_Y' 'Pos_Z' 'Forward_X' 'Forward_Y' 'Forward_Z'});

dir_23   = [(Cone_pos(3,1) - Cone_pos(2,1)) (Cone_pos(3,3) - Cone_pos(2,3))];
dir_trig = [(Trig_pos(1,1) - Cone_pos(3,1)) (Trig_pos(1,3) - Cone_pos(3,3))];
f_return_angle = anglebetween(dir_23,dir_trig);

if(~isfield(Group,'Reconstructed'))
    disp('Please run CalculateTrckingPath first');
    return;
end

cond = Group.CondTable{1,pId}.Condition(trialId);
condStr = '';

if(cond == 1)
    condStr = 'No change';
elseif(cond == 2)
    condStr = 'No distal cues';
else
    condStr = 'No optic flow';
end

calculatedAngle = Group.Reconstructed{1,pId}.InboundBodyRotation(trialId);
real_return_angle = Group.Reconstructed{1,pId}.RealReturnAngle(trialId);

close all; clc;

% Extracting after reaching cone 3
if(config.cutFromConeThree)
    Extracted_pos = Extracted_pos(Extracted_pos.Time > Group.FlagTrigTimes{1,pId}{trialId,3} - 2,:);
end

% Plotting according to the speed
CreateCustomFigure;

hold all

ax = gca;
ax.FontSize = 20;
xlim([-1 4])
ylim([-1 4])

xlabel('x (m)')
ylabel('y (m)')

title([...
    'Participant: ', num2str(pId),...
    ' Trial: ', num2str(trialId),...
    ' Body Rotation: ',num2str(calculatedAngle,'%.0f') ,...
    ' Inferred Angle: ', num2str(f_return_angle,'%.0f'),...
    ' Real Inferred Angle: ' , num2str(real_return_angle,'%.0f')...
    ' Condition: ' , condStr]...
    ,FontSize=20);


axis square

tracking_size = height(Extracted_pos);

% Plot the position of the third cone
pCone1 = plot(Cone_pos(1,1),Cone_pos(1,3),'Marker','d','LineStyle','none','MarkerFaceColor',config.color_scheme_npg(2,:),'MarkerEdgeColor','black','MarkerSize',18); 
pCone2 = plot(Cone_pos(2,1),Cone_pos(2,3),'Marker','d','LineStyle','none','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',18); 
pCone3 = plot(Cone_pos(3,1),Cone_pos(3,3),'Marker','d','LineStyle','none','MarkerFaceColor',config.color_scheme_npg(1,:),'MarkerEdgeColor','black','MarkerSize',18);
pTrigPos = plot(Trig_pos(1,1),Trig_pos(1,3),'Marker','x','LineStyle','none','MarkerFaceColor',config.color_scheme_npg(2,:),'MarkerEdgeColor','black','MarkerSize',18);

leg = legend([pCone1 pCone2 pCone3 pTrigPos], 'Cone 1', 'Cone 2', 'Cone 3', 'Response','AutoUpdate','off');
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

    if(k > 2 && k < tracking_size)
        delete(qv);
    end

end

hold off

end