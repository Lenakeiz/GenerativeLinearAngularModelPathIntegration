function VisualizeRealtimeTrackingData(Group, pId, trialId, playbackSpeed, varargin)

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

if(~isfield(Group,'Reconstructed'))
    disp('Please run CalculateTrckingPath first');
    return;
end

cond = Group.CondTable{1,pId}.Condition(trialId);
outofbound = Group.CondTable{1,pId}.OutOfBound(trialId);

condStr = '';
if(cond == 1)
    condStr = 'no change';
elseif(cond == 2)
    condStr = 'no distal cues';
else
    condStr = 'no optic flow';
end

outStr = '';
if(outofbound == 1)
    outStr = 'yes';
else
    outStr = 'no';
end

calculatedAngle = Group.Reconstructed{1,pId}.InboundBodyRotation(trialId);
real_return_angle = Group.Reconstructed{1,pId}.RealReturnAngle(trialId);
f_return_angle = Group.Reconstructed{1,pId}.InferredReturnAngle(trialId);

close all; clc;

OoB = [100 0 100];
ReconstrutedOoB = [100 0 100];
if (outofbound == 1)
    OoB = Group.OutOfBoundPos{1,pId}{trialId,1};
    if(isfield(Group,'ReconstructedOOB'))
        ReconstrutedOoB = Group.ReconstructedOOB{1,pId}.ReconstructedOoB{trialId,1};
    end    
end

% Extracting after reaching cone 3
if(config.cutFromConeThree)
    Extracted_pos = Extracted_pos(Extracted_pos.Time > Group.FlagTrigTimes{1,pId}{trialId,3} - 2,:);
end

% Plotting according to the speed
CreateCustomFigure;

hold all

ax = gca;
ax.FontSize = 20;
ax.FontName = 'Times New Roman';

% Dinamically calculate bounds
maxX = max([Cone_pos(:,1); Trig_pos(1,1)]);
maxY = max([Cone_pos(:,3); Trig_pos(1,3)]);
minX = min([Cone_pos(:,1); Trig_pos(1,1)]);
minY = min([Cone_pos(:,3); Trig_pos(1,3)]);
absMax = ceil(max(maxX, maxY));
absMin = floor(min(minX,minY));
offSetFromMax = 0.5;

xlim([absMin-offSetFromMax absMax+offSetFromMax]);
ylim([absMin-offSetFromMax absMax+offSetFromMax]);

xlabel('x (m)')
ylabel('y (m)')

title(['Path visualizer'],FontSize=20);

axis square

boxwidth = [0.15 0.2];
leftCorner = [0.7 0.4];
dim = [leftCorner boxwidth];
str = {...
    ['\textbf{Trial Info}'],...
    ['\textit{Participant:} ', num2str(pId) ],...
    ['\textit{Trial:} ', num2str(trialId)],...
    ['\textit{Body Rotation:} $',num2str(calculatedAngle,'%.0f'),'^{\circ}$'],...
    ['\textit{Inferred Angle:} $', num2str(f_return_angle,'%.0f'),'^{\circ}$'],...
    ['\textit{Real Inferred Angle:} $' , num2str(real_return_angle,'%.0f'),'^{\circ}$'],...
    ['\textit{Condition:} ' , condStr],...
    ['\textit{Out of Bound:} ', outStr]...
    };

annotation('textbox',dim,String=str,FontSize=15,Interpreter='latex', FitBoxToText='on');

tracking_size = height(Extracted_pos);

% Plot the position of the third cone
pCone1 = plot(Cone_pos(1,1),Cone_pos(1,3),'Marker','d','LineStyle','none','MarkerFaceColor',config.color_scheme_npg(2,:),'MarkerEdgeColor','black','MarkerSize',18); 
pCone2 = plot(Cone_pos(2,1),Cone_pos(2,3),'Marker','d','LineStyle','none','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','black','MarkerSize',18); 
pCone3 = plot(Cone_pos(3,1),Cone_pos(3,3),'Marker','d','LineStyle','none','MarkerFaceColor',config.color_scheme_npg(1,:),'MarkerEdgeColor','black','MarkerSize',18);
pOoB = plot(OoB(1,1),OoB(1,3),'Marker','*',Color=config.color_scheme_npg(4,:),MarkerSize=18, LineWidth=3);
pRecconstructedOoB = plot(ReconstrutedOoB(1,1),ReconstrutedOoB(1,3),'Marker','*',Color=config.color_scheme_npg(7,:),MarkerSize=18,LineWidth=3);
pTrigPos = plot(Trig_pos(1,1),Trig_pos(1,3),'Marker','x','LineStyle','none','MarkerFaceColor',config.color_scheme_npg(2,:),'MarkerEdgeColor','black','MarkerSize',18,LineWidth=3);

leg = legend([pCone1 pCone2 pCone3 pTrigPos pOoB pRecconstructedOoB], 'Cone 1', 'Cone 2', 'Cone 3','Response','Out of Bound(OoB)','Reconstructed OoB','AutoUpdate','off');
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