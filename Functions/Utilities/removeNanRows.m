function [NonNanDat, nanIdx] = removeNanRows(Dat)
%remove Nan Rows in Dat (2D array)
nanIdx = isnan(sum(Dat,2));
NonNanDat = Dat(~nanIdx,:);
