function NonNanDat = removeNanRows(Dat)
%remove Nan Rows in Dat (2D array)
nanIdx = isnan(Dat(:,1));
NonNanDat = Dat(~nanIdx,:);
