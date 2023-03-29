function [Y, GroupNames, ConditionNames]=ReGroupData(Data,groupname)
%% ReGroupData
% Andrea Castegnaro, UCL, uceeaca@ucl.ac.uk
% Group all the group data from three environmental conditions into a long
% numeric vector for further two-way anova analysis and also output the
% factor names
% ===================================================================================

% Removing nans
isnan_idx = ~isnan(Data(1,:));
nochange = Data(1,isnan_idx);
isnan_idx = ~isnan(Data(2,:));
nodistalcues = Data(1,isnan_idx);
isnan_idx = ~isnan(Data(3,:));
noopticflow = Data(1,isnan_idx);

Y = [nochange, nodistalcues, noopticflow];

numSubjsNoChange     = size(nochange,2);
numSubjsNoDistalCues = size(nodistalcues,2);
numSubjsNoOpticFlow  = size(noopticflow,2);

GroupNames = [string(repmat({groupname},1,numSubjsNoChange)),...
              string(repmat({groupname},1,numSubjsNoDistalCues)),...
              string(repmat({groupname},1,numSubjsNoOpticFlow))];

ConditionNames = [string(repmat({'NoChange'},1 ,numSubjsNoChange)),...
                  string(repmat({'NoDistalCue'},1,numSubjsNoDistalCues)),...
                  string(repmat({'NoOpticFlow'},1,numSubjsNoOpticFlow))];  

end