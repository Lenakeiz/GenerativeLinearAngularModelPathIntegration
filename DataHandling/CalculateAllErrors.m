function [ RetData ] = CalculateAllErrors(GroupData)
%% CalculateAllErrors
% Function wrapper that will calculate the absolute distance error on the
% inbound path for each participant of a given group. See nested
% CalculateErrors function for details on the calculation.
% ===================================================================================

sampleSize = size(GroupData.FlagPos,2);

for t = 1:sampleSize
    [Errors{t} OoBInfo{t}] = CalculateErrors(GroupData.FlagPos{t}, GroupData.TrigPos{t}, GroupData.OutOfBoundPos{t}, GroupData.CondTable{t});
end

RetData.Errors = Errors;
RetData.OoBInfo = OoBInfo;

end

