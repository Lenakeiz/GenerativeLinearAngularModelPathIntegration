function VisualizeRealtimeTrackingData(Group, pId, trialId, playbackSpeed, varargin)
%% VisualizeRealtimeTrackingData
% Andrea Castegnaro, UCL, 2022 uceeaca@ucl.ac.uk
% Function used to manuyally visualize the reconstructed data from a
% participant.
% Group: group from which we are going to visualize the trial
% pId: identify participant within the group
% trialId: identify trial within participant
% playbackSpeed: set the speed for the video
% varargin: 'cutconethree' - extract data from cone 3 only, 'videoplot' -
% draw the plot as video or not, 'cutsecbeforecone1' - number of seconds to
% track before reaching cone 1, 'displaytrialinfo' - plot a legend with 
% This function has been used to create Figure 1b.
% ===================================================================================

% Anonymous function to calculate angle between two vectors.
anglebetween = @(v,w) atan2d(w(:,2).*v(:,1) - v(:,2).*w(:,1), v(:,1).*w(:,1) + v(:,2).*w(:,2));

% Color pattern for our paper
ColorPattern;

% Preparing output
config.ResultFolder = pwd + "/Output/Fig1b";

if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

% Creating figure
figureOpen = figure('visible','off','Position', [100 100 850 500]);
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultTextFontSize',16)  

% Reconstructing a visualization of participants walk and head drection after reaching the third
% cone
config.cutFromConeThree = false;
% Animate animation or not
config.videoplot = true;
% Seconds to start tracking before reaching cone 1 (-100 means get all
% information in tracking path from the current trial. Current trial starts
% when the participants ended the previous one.
config.cutseconds = -100; 
% Diplay a legend with all the trial info. Necessary for the manual
% checking of the data
config.displayTrialInfo = false; 

if exist('varargin','var')
    for i = 1:2:nargin-4
        if(strcmpi(varargin{i},'cutconethree'))
            config.cutFromConeThree = varargin{i+1};
        end
        if(strcmpi(varargin{i},'videoplot'))
            config.videoplot = varargin{i+1};
        end
        if(strcmpi(varargin{i},'cutsecbeforecone1'))
            config.cutseconds = varargin{i+1};
        end
        if(strcmpi(varargin{i},'displaytrialinfo'))
            config.displayTrialInfo = varargin{i+1};
        end
    end
end

Cone_pos      = Group.FlagPos{1,pId}{trialId,1};
Trig_pos      = Group.TrigPos{1,pId}{trialId,1};
Extracted_pos = Group.Path{1,pId}{trialId,1};
Extracted_pos = array2table(Extracted_pos,"VariableNames",{'Time' 'Pos_X' 'Pos_Y' 'Pos_Z' 'Forward_X' 'Forward_Y' 'Forward_Z'});

if(~isfield(Group,'Reconstructed'))
    disp('Please run CalculateTrackingPath function on the group before this one');
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

% Information produced in CalculateTrackingPath
calculatedAngle = Group.Reconstructed{1,pId}.InboundBodyRotation(trialId);
real_return_angle = Group.Reconstructed{1,pId}.RealReturnAngle(trialId);
f_return_angle = Group.Reconstructed{1,pId}.InferredReturnAngle(trialId);

% Adding information if this is an out of bound trial
OoB = [nan 0 nan];
ReconstrutedOoB = [nan 0 nan];
if (outofbound == 1)
    OoB = Group.OutOfBoundPos{1,pId}{trialId,1};
    if(isfield(Group,'ReconstructedOOB'))
        ReconstrutedOoB = Group.ReconstructedOOB{1,pId}.ReconstructedOoB{trialId,1};
    end    
end

% Extracting tracking information after cone three has been reached
if(config.cutFromConeThree)
    Extracted_pos = Extracted_pos(Extracted_pos.Time > Group.FlagTrigTimes{1,pId}{trialId,3} - 2,:);
end
if(config.cutseconds > 0)
    Extracted_pos = Extracted_pos(Extracted_pos.Time > Group.FlagTrigTimes{1,pId}{trialId,1} - config.cutseconds,:);
end

hold on

ax = gca;
ax.FontSize = 20;
ax.FontName = 'Arial';

% Dinamically calculate bounds for the plot
maxX = max([Cone_pos(:,1); Trig_pos(1,1); OoB(1,1)]);
maxY = max([Cone_pos(:,3); Trig_pos(1,3); OoB(1,3)]);
minX = min([Cone_pos(:,1); Trig_pos(1,1); OoB(1,1)]);
minY = min([Cone_pos(:,3); Trig_pos(1,3); OoB(1,3)]);
absMax = ceil(max(maxX, maxY));
absMin = floor(min(minX,minY));
offSetFromMax = 0.5;

xlabel('x (m)')
ylabel('y (m)')

% Creating a legend with information calculated in CalculateTrackingPath
boxwidth = [0.15 0.2];
leftCorner = [0.65 0.4];
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

if(config.displayTrialInfo)
    ann = annotation('textbox',dim,String=str,FontSize=15,Interpreter='latex', FitBoxToText='on', FontUnits='points');
    ann.FontName = "Arial";
    ann.FontSize = 14;
end

tracking_size = height(Extracted_pos);

% Plot the position of the cones
pCone1 = plot(Cone_pos(1,1),Cone_pos(1,3),'Marker','d','LineStyle','none','MarkerFaceColor',config.color_scheme_npg(2,:),'MarkerEdgeColor','black','MarkerSize',18); 
pCone2 = plot(Cone_pos(2,1),Cone_pos(2,3),'Marker','d','LineStyle','none','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','black','MarkerSize',18); 
pCone3 = plot(Cone_pos(3,1),Cone_pos(3,3),'Marker','d','LineStyle','none','MarkerFaceColor',config.color_scheme_npg(1,:),'MarkerEdgeColor','black','MarkerSize',18);
% Mark out of bound position when available
pOoB = plot(OoB(1,1),OoB(1,3),'Marker','*',Color=config.color_scheme_npg(4,:),MarkerSize=18, LineWidth=3);
% Mark the reconstructed out of bound position (this position is set to a
% distance for which it will not weight on the mean)
pRecconstructedOoB = plot(ReconstrutedOoB(1,1),ReconstrutedOoB(1,3),'Marker','*',Color=config.color_scheme_npg(7,:),MarkerSize=18,LineWidth=3);
% Plot the triggered position
pTrigPos = plot(Trig_pos(1,1),Trig_pos(1,3),'Marker','x','LineStyle','none','MarkerFaceColor',config.color_scheme_npg(2,:),'MarkerEdgeColor','black','MarkerSize',18,LineWidth=3);

if(config.videoplot)
    % Draw frame by frame
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
        
        drawnow; pause(playbackSpeed);
    
        if(k > 2 && k < tracking_size)
            delete(qv);
        end        
    end
else
    tempDir = [Extracted_pos.Forward_X(2) Extracted_pos.Forward_Z(2)];
    tempDir = tempDir/norm(tempDir);
    % Making participant orientation smaller just for visualization
    tempDir = tempDir.*0.5;

    qv = quiver(Extracted_pos.Pos_X(2),Extracted_pos.Pos_Z(2),tempDir(1), tempDir(2), 'off');
    qv.Color = config.color_scheme_npg(5,:);
    qv.LineWidth = 5;
    qv.MaxHeadSize = 2;

    tempDir = [Extracted_pos.Forward_X(tracking_size) Extracted_pos.Forward_Z(tracking_size)];
    tempDir = tempDir/norm(tempDir);
    tempDir = tempDir.*0.5;

    qv = quiver(Extracted_pos.Pos_X(tracking_size),Extracted_pos.Pos_Z(tracking_size),tempDir(1), tempDir(2), 'off');
    qv.Color = config.color_scheme_npg(5,:);
    qv.LineWidth = 5;
    qv.MaxHeadSize = 2;

    trackDataPlot = plot(Extracted_pos.Pos_X,Extracted_pos.Pos_Z,'Marker','none','LineStyle','-','LineWidth',2.5,'Color',config.color_scheme_npg(4,:));
end

if(outofbound == 1)
    leg = legend([pCone1 pCone2 pCone3 pTrigPos pOoB pRecconstructedOoB], 'Cone 1', 'Cone 2', 'Cone 3','Response','Out of Bound(OoB)','Reconstructed OoB','AutoUpdate','off');
else
    leg = legend([pCone1 pCone2 pCone3 pTrigPos trackDataPlot], 'Cone 1', 'Cone 2', 'Cone 3','Response','Tracking data','AutoUpdate','off');
end
leg.Location = 'northeast';
leg.FontSize = 15;

axis square

xlim([absMin-offSetFromMax absMax+offSetFromMax]);
ylim([absMin-offSetFromMax absMax+offSetFromMax]);
yticks(xticks);

ax = gca;
ax.LineWidth = 2.0;

hold off

% Export figure
exportgraphics(figureOpen,config.ResultFolder+"/ExampleTrackingPath.png",'Resolution',300);
exportgraphics(figureOpen,config.ResultFolder+"/ExampleTrackingPath.pdf",'Resolution',300, 'ContentType','vector');
disp("Save complete!");
end