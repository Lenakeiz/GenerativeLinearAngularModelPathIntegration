function [MCI, StatusIndex]= MergeMCI(MCIPos, MCINeg, MCIUnk)
%% MergeMCI
% Merging all of the MCI data into a single structure. Useful for running
% ANOVAS and grouped statistics.
% Output:
% MCI - merged structure
% Status Index - MCI identifier coded as follows: positive = 2, negative = 1, unknown = 0
% ===================================================================================

    numConds = size(MCIPos,2);
    StatusIndex = 2*ones(size(MCIPos{1},1),1);
    StatusIndex = [StatusIndex; ones(size(MCINeg{1},1),1)];
    StatusIndex = [StatusIndex; zeros(size(MCIUnk{1},1),1)];
    MCI = cell(1, numConds);
    for condidx =1:numConds
        MCI{condidx} = [MCIPos{condidx};MCINeg{condidx};MCIUnk{condidx}];
    end
end