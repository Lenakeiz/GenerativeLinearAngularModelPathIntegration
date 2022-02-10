%% TransformPaths.m
% Transform the paths so they are all facing the real axis and starting from (0,0).
% Also mirrors all the left turns in the outbound path to have all paths
% being rotated towards the right
% 
function [outData] = TransformPaths(Data)
%TRANSFORMPATHS Summary of this function goes here
%   Detailed explanation goes here
outData = Data;

% SET to debug plot
PLOT_FEEDBACK = 0;

anglebetween = @(va,vb) atan2d(va(:,1).*vb(:,2) - va(:,2).*vb(:,1), va(:,1).*vb(:,1) + va(:,2).*vb(:,2));

sampleSize = size(Data.FlagPos,2);

if(PLOT_FEEDBACK == 1)
close all;
figure;
end

for npart = 1:sampleSize
        
    flagpos = Data.FlagPos{npart};
    trigpos = Data.TrigPos{npart};
    
    flagposTemp = cell(0);
    trigposTemp = cell(0);
    
    %Looping through the trials
    for tr = 1:length(flagpos)
        X_coord = [flagpos{tr}(:,[1,3]);trigpos{tr}([1,3])];
        % Translating to zero
        X_coord_zero = X_coord - X_coord(1,:);
        
        realAxisDirection = [0 0; 1 0];
        angle = anglebetween(realAxisDirection,X_coord_zero(1:2,:));
        angle = angle(2);
        
        % Rotation matrix in 2D
        RMatrix = [cosd(angle) sind(angle); -sind(angle) cosd(angle)];
        
        % Applying rotation matrix
        X_coord_rotated = (RMatrix * X_coord_zero')';
        
        % Calculating internal angle between the first and second segment
        firstSegment  = X_coord_rotated(2,:) - X_coord_rotated(1,:);
        secondSegment = X_coord_rotated(3,:) - X_coord_rotated(2,:);
        internalAngle = anglebetween(-firstSegment,secondSegment);
        
        % Mirroring left turns means inverting the y-position of segments
        % after the second one
        if (internalAngle > 0)
            X_coord_rotated(3,2) = X_coord_rotated(3,2) * -1;
            X_coord_rotated(4,2) = X_coord_rotated(4,2) * -1;
        end
        
        %internalAngle = anglebetween(X_coord_rotated(1:2,:),X_coord_zero(2:3,:));
        %internalAngle = internalAngle(2);
        
        %Adding back a dummy y coordinate which is the height of the
        %participant and should never be used
        flagposTemp{tr,1} = [X_coord_rotated(1:3,1) ones(3,1)*2 X_coord_rotated(1:3,2)];
        trigposTemp{tr,1} = [X_coord_rotated(4,1) ones(1,1)*2 X_coord_rotated(4,2)];

        if(PLOT_FEEDBACK == 1)
            subplot(1,2,1);
            ppData = plot(X_coord_zero(:,1),X_coord_zero(:,2),'k'); hold on;
            ppData.LineWidth =  1.5;

            ttex = text(X_coord_zero(:,1),X_coord_zero(:,2),['0' '1' '2' '3']','FontSize',20);
                        
            ttex = text(2.5,2.5,num2str(internalAngle),'FontSize',20);
            ttex.FontSize = 20;

            xlim([-3 3]);
            ylim([-3 3]);
            axis square;
            hold off;
            drawnow;

            subplot(1,2,2);
            ppData = plot(X_coord_rotated(:,1),X_coord_rotated(:,2),'k'); hold on;
            ppData.LineWidth =  1.5;
            
            % Plotting indices of the transformed vertices
            ttex = text(X_coord_rotated(:,1),X_coord_rotated(:,2),['0' '1' '2' '3']','FontSize',20);
                        
            ttex = text(2.5,2.5,num2str(angle),'FontSize',20);
            ttex.FontSize = 20;

            xlim([-3 3]);
            ylim([-3 3]);
            axis square;

            hold off;
            drawnow;

            waitforbuttonpress;
        end
        
    end
    
    outData.FlagPos{npart} = flagposTemp;
    outData.TrigPos{npart} = trigposTemp;
    
end

end

