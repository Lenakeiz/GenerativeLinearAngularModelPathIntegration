function [Y, GroupNames, ConditionNames]=ReGroupData(Data,groupname)
%% ReGroupData
% Andrea Castegnaro, UCL, uceeaca@ucl.ac.uk
% Group all the group data from three environmental conditions into a long
% numeric vector for further two-way anova analysis and also output the
% factor names
% ===================================================================================

Y = [Data(1,:), Data(2,:), Data(3,:)];

numSubjs = size(Data,2);

GroupNames = [string(repmat({groupname},1,numSubjs)),...
              string(repmat({groupname},1,numSubjs)),...
              string(repmat({groupname},1,numSubjs))];

ConditionNames = [string(repmat({'NoChange'},1,numSubjs)),...
                  string(repmat({'NoDistalCue'},1,numSubjs)),...
                  string(repmat({'NoOpticFlow'},1,numSubjs))];  

end