function [outData] = TransformPaths(Data)
%% TransformPaths
% Andrea Castegnaro, UCL, 2022 uceeaca@ucl.ac.uk
% Transform all the paths so they are all facing the positive x-axis and
% starting from (0,0). Mirrors all the right turns in the outbound
% so that all paths are left turns.
% Once flipped calculate the errors for each of trial.
% Input structure "Data" is group data loaded from the .mat file in Data
% folder and is being overriden as an output
% ===================================================================================

% Overriding the outout data
outData = Data;

% In function debug plot. Set to one to plot output.
PLOT_FEEDBACK = 0;

% Anonymous function to calculate angle between two vectors.
anglebetween = @(va,vb) atan2d(va(:,1).*vb(:,2) - va(:,2).*vb(:,1), va(:,1).*vb(:,1) + va(:,2).*vb(:,2));

sampleSize = size(Data.FlagPos,2);

if(PLOT_FEEDBACK == 1)
    close all;
    figure;
end

for npart = 1:sampleSize

    flagpos       = Data.FlagPos{npart};
    trigpos       = Data.TrigPos{npart};
    outofboundpos = Data.OutOfBoundPos{npart};

    flagposTemp       = cell(0);
    trigposTemp       = cell(0);
    outofboundpostemp = cell(0);

    if(isempty(outofboundpos{1}))
        % This participant will be subsequenlty discarded from the model
        % fitting as we could not reconstruct the angles in the OoB trials
        disp("Participant #" + npart + " does not have a registered OoB position");
    end

    track_flipping = [];

    for tr = 1:length(flagpos)

        flip = 0;
        firstConePositions  = flagpos{tr}(1 : 3 : end, [1,3]);
        secondConePositions = flagpos{tr}(2 : 3 : end, [1,3]);
        thirdConePositions  = flagpos{tr}(3 : 3 : end, [1,3]);
        
        % Tracking data
        path = Data.Path{1,npart}{tr,1};
        % Triangle coordinates
        X_coord = [flagpos{tr}(:,[1,3]);trigpos{tr}([1,3])];
        
        % Step 1: Translate to zero coordinate
        X_coord_zero = X_coord - X_coord(1,:);

        if(~isempty(outofboundpos{tr}))
            X_outOfBound = outofboundpos{tr}([1,3]);
        else
            X_outOfBound = nan(1,2);
        end

        X_outOfBound_zero = X_outOfBound - X_coord(1,:);
        path(:,[2 4]) = path(:,[2 4]) - X_coord(1,:);

        % Step 2:  Calculate angle between positive x-axis and l1 segment.
        % Apply rotation matrix to align data with positive x-axis.
        realAxisDirection = [0 0; 1 0];
        angle = anglebetween(realAxisDirection,X_coord_zero(1:2,:));
        angle = angle(2);
        RMatrix = [cosd(angle) sind(angle); -sind(angle) cosd(angle)];

        X_coord_rotated = (RMatrix * X_coord_zero')';
        X_outOfBound_rotated = (RMatrix * X_outOfBound_zero')';
        

        path(:,[2 4]) = (RMatrix * path(:,[2 4])')';
        path(:,[5 7]) = (RMatrix * path(:,[5 7])')';

        % Step 3: flip right turn trials
        firstSegment  = X_coord_rotated(2,:) - X_coord_rotated(1,:);
        secondSegment = X_coord_rotated(3,:) - X_coord_rotated(2,:);
        internalAngle = anglebetween(-firstSegment,secondSegment);

        if (internalAngle > 0)
            flip = 1;
            X_coord_rotated(3,2) = X_coord_rotated(3,2) * -1;
            X_coord_rotated(4,2) = X_coord_rotated(4,2) * -1;
            X_outOfBound_rotated(1,2) = X_outOfBound_rotated(1,2) * -1;
            path(:,[4]) = path(:,[4]).*-1;
            path(:,[7]) = path(:,[7]).*-1;
        end

        % Marking flipped trials
        track_flipping = [track_flipping;flip];
        
        % For consistency with our data format we reinsert back the Y
        % coordinate (which is the height of cones/ headset). We use a dummy variable 
        flagposTemp{tr,1} = [X_coord_rotated(1:3,1) ones(3,1)*2 X_coord_rotated(1:3,2)];
        trigposTemp{tr,1} = [X_coord_rotated(4,1) ones(1,1)*2 X_coord_rotated(4,2)];
        outofboundpostemp{tr,1} = [X_outOfBound_rotated(1,1) 1.0 X_outOfBound_rotated(1,2)];

        % Save to the output variable
        outData.Path{1,npart}{tr,1} = path;

        % In function quick visualization 
        if(PLOT_FEEDBACK == 1 && Data.CondTable{npart}.OutOfBound(tr) == 1)

            disp("Participant " + npart + " trial " +  tr);
            disp("Out of Bound Pos: " + X_outOfBound);
            disp("Out of Bound Rec Pos: " + X_outOfBound_rotated);

            %Plotting real data
            subplot(2,2,1);
            a = [firstConePositions;secondConePositions;thirdConePositions;trigpos{tr}([1,3]);X_outOfBound];
            ppData = plot(X_coord(:,1),X_coord(:,2),'k'); hold on;
            ppData = plot(a(:,1),a(:,2)); ppData.LineWidth =  1.2;
            ttex = text(a(:,1),a(:,2),["0" "1" "2" "3" "OoB"]','FontSize',15);
            title("Real Data (OoB)");
            xlim([-4 4]); ylim([-4 4]); axis square; hold off;

            %Plotting data centred on zero
            subplot(2,2,2);
            a = [X_coord_zero; X_outOfBound_zero];
            ppData = plot(a(:,1),a(:,2),'k'); hold on;
            ppData.LineWidth =  1.2; ttex = text(a(:,1),a(:,2),["0" "1" "2" "3" "OoB"]','FontSize',15);
            title("Real Data Zero Cent(OoB)");
            xlim([-4 4]); ylim([-4 4]); axis square; hold off;

            %Plotting rotated and transformed data
            subplot(2,2,3);
            a=[X_coord_rotated;X_outOfBound_rotated];
            ppData = plot(a(:,1),a(:,2),'k'); hold on;
            ppData.LineWidth =  1.2;
            ttex = text(a(:,1),a(:,2),["0" "1" "2" "3" "OoB"]','FontSize',15);
            title("Real Data Zero Cent(OoB)");
            xlim([-4 4]); ylim([-4 4]); axis square; hold off; 

            drawnow; waitforbuttonpress;

        end
    end

    % Save information on the output structure
    outData.Reconstructed{1,npart}.FlippedTrial = track_flipping;

    outData.FlagPos{npart} = flagposTemp;
    outData.TrigPos{npart} = trigposTemp;
    outData.OutOfBoundPos{npart} = outofboundpostemp;

end

ErrorInfo = CalculateAllErrors(outData);
outData.Errors = ErrorInfo.Errors;
outData.ReconstructedOOB = ErrorInfo.OoBInfo;

end