function MCI = MergeMCI(MCIPos, MCINeg, MCIUnk)
    numConds = size(MCIPos,2);
    MCI = cell(1, numConds);
    for condidx =1:numConds
        MCI{condidx} = [MCIPos{condidx};MCINeg{condidx};MCIUnk{condidx}];
    end
end