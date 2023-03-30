function [ RetData ] = CalculateAllOoBTrials(GroupData)
%% CalculateAllOoBTrials
% Andrea Castegnaro, UCL, 2022 uceeaca@ucl.ac.uk
% Function wrapper that will preprocess the Out of Bound trials for a given
% group. See nested CalculateOoB function for details and online methods.
% ===================================================================================

sampleSize = size(GroupData.FlagPos,2);

for t = 1:sampleSize
    [OoBInfo{t}] = CalculateOoB(GroupData.FlagPos{t}, GroupData.TrigPos{t}, GroupData.OutOfBoundPos{t}, GroupData.CondTable{t});
end

RetData.OoBInfo = OoBInfo;

end

