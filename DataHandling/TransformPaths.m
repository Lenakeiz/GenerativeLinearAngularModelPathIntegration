function [outData] = TransformPaths(Data)
% TRANSFORMPATHS Transform all the paths so they are all facing the real
% axis and starting from (0,0). Mirrors all the left turns in the outbound
% path to have all paths being right turn. In that case it s easier to
% understand what s the overshooting if using the formula for angles in
% anglebetween. 
% Input structure "Data" is being overriden and set as an output
% ===================================================================================
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
        
    flagpos       = Data.FlagPos{npart};
    trigpos       = Data.TrigPos{npart};
    outofboundpos = Data.OutOfBoundPos{npart};
    
    flagposTemp       = cell(0);
    trigposTemp       = cell(0);
    outofboundpostemp = cell(0);

    if(isempty(outofboundpos{1}))
        disp("Participant #" + npart + " does not have a registered OuB position");
    end

    track_flipping = [];

    %Looping through the trials
    for tr = 1:length(flagpos)

        flip = 0;
        firstConePositions  = flagpos{tr}(1 : 3 : end, [1,3]);
        secondConePositions = flagpos{tr}(2 : 3 : end, [1,3]);
        thirdConePositions  = flagpos{tr}(3 : 3 : end, [1,3]);
        
        X_coord = [flagpos{tr}(:,[1,3]);trigpos{tr}([1,3])];
        % Translating to zero
        X_coord_zero = X_coord - X_coord(1,:);
        
        if(~isempty(outofboundpos{tr}))
            X_outOfBound = outofboundpos{tr}([1,3]);
        else            
            X_outOfBound = nan(1,2);
        end

        X_outOfBound_zero = X_outOfBound - X_coord(1,:);

        realAxisDirection = [0 0; 1 0];
        angle = anglebetween(realAxisDirection,X_coord_zero(1:2,:));
        angle = angle(2);
        
        % Rotation matrix in 2D
        RMatrix = [cosd(angle) sind(angle); -sind(angle) cosd(angle)];
        
        % Applying rotation matrix
        X_coord_rotated = (RMatrix * X_coord_zero')';
        X_outOfBound_rotated = (RMatrix * X_outOfBound_zero')';
        % Calculating internal angle between the first and second segment
        firstSegment  = X_coord_rotated(2,:) - X_coord_rotated(1,:);
        secondSegment = X_coord_rotated(3,:) - X_coord_rotated(2,:);
        internalAngle = anglebetween(-firstSegment,secondSegment);
        
        % Mirroring left turns means inverting the y-position of segments
        % after the second one
        if (internalAngle > 0)
            flip = 1;
            X_coord_rotated(3,2) = X_coord_rotated(3,2) * -1;
            X_coord_rotated(4,2) = X_coord_rotated(4,2) * -1;
            X_outOfBound_rotated(1,2) = X_outOfBound_rotated(1,2) * -1;
        end
        
        track_flipping = [track_flipping;flip];
        %internalAngle = anglebetween(X_coord_rotated(1:2,:),X_coord_zero(2:3,:));
        %internalAngle = internalAngle(2);
        
        %Adding back a dummy y coordinate which is the height of the
        %participant and should never be used
        flagposTemp{tr,1} = [X_coord_rotated(1:3,1) ones(3,1)*2 X_coord_rotated(1:3,2)];
        trigposTemp{tr,1} = [X_coord_rotated(4,1) ones(1,1)*2 X_coord_rotated(4,2)];
        outofboundpostemp{tr,1} = [X_outOfBound_rotated(1,1) 1.0 X_outOfBound_rotated(1,2)];

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
            %ttex = text(2.5,2.5,num2str(internalAngle),'FontSize',20);
            %ttex.FontSize = 20; 
            xlim([-4 4]); ylim([-4 4]); axis square; hold off; drawnow; waitforbuttonpress;

            %Plotting data centred on zero
            subplot(2,2,2);
            a = [X_coord_zero; X_outOfBound_zero];
            ppData = plot(a(:,1),a(:,2),'k'); hold on; 
            ppData.LineWidth =  1.2; ttex = text(a(:,1),a(:,2),["0" "1" "2" "3" "OoB"]','FontSize',15);
            title("Real Data Zero Cent(OoB)");
%           ttex = text(2.5,2.5,num2str(angle),'FontSize',20);
%           ttex.FontSize = 20;
            xlim([-4 4]); ylim([-4 4]); axis square; hold off; drawnow; waitforbuttonpress;

            %Plotting rotated and transformed data
            subplot(2,2,3);
            a=[X_coord_rotated;X_outOfBound_rotated];
            ppData = plot(a(:,1),a(:,2),'k'); hold on; 
            ppData.LineWidth =  1.2;
            ttex = text(a(:,1),a(:,2),["0" "1" "2" "3" "OoB"]','FontSize',15);
            title("Real Data Zero Cent(OoB)");
%           ttex = text(2.5,2.5,num2str(angle),'FontSize',20);
%           ttex.FontSize = 20;
            xlim([-4 4]); ylim([-4 4]); axis square; hold off; drawnow; waitforbuttonpress;
        end        
    end
    
    outData.Reconstructed{1,npart}.FlippedTrial = track_flipping;

    outData.FlagPos{npart} = flagposTemp;
    outData.TrigPos{npart} = trigposTemp;
    outData.OutOfBoundPos{npart} = outofboundpostemp;

end

ErrorInfo = CalculateAllErrors(outData);
outData.Errors = ErrorInfo.Errors;
outData.ReconstructedOOB = ErrorInfo.OoBInfo;

end