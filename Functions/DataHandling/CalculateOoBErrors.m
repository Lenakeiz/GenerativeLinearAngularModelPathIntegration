function [ RetData ] = CalculateOoBErrors(GroupData)
    % This script will calculate a set of errors to assess whether the 
    % errors
    sampleSize = size(GroupData.FlagPos,2);

    for t = 1:sampleSize   
        %Getting abs error distance
        [Errors{t} OoBInfo{t}] = CalculateOoB(GroupData.FlagPos{t}, GroupData.TrigPos{t}, GroupData.OutOfBoundPos{t}, GroupData.CondTable{t});
    end

    RetData.Errors = Errors;
    RetData.OoBInfo = OoBInfo;
    
end


