function [NonNanDat, nanIdx] = removeNanRows(Dat)
%% removeNanRows
% Andrea Castegnaro, UCL, uceeaca@ucl.ac.uk
% Remove any nan rows from Dat array
% ===================================================================================
nanIdx = isnan(sum(Dat,2));
NonNanDat = Dat(~nanIdx,:);
