% Get pixel position of monitors.
close all;

fh = figure('units', 'pixels');
MP = get(0, 'MonitorPositions');
N = size(MP, 1);
newPosition = MP(1,:);
if size(MP, 1) == 1

else
    newPosition(1) = newPosition(1) + MP(N,1);
end
fh.set('Position', newPosition, 'units', 'normalized');
fh.WindowState = 'maximized';

set(gcf,'color','w');

clear MP N newPosition
%