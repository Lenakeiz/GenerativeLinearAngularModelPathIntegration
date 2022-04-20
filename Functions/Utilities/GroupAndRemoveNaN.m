function [Y, GroupNames, ConditionNames]=GroupAndRemoveNaN(Data,param_idx, groupname)
%Group all the data from all groups and all conditions into a long
%numeric vector for further two-way anova analysis 
%and also output the names 

Y = [Data{1}(:,param_idx)',Data{2}(:,param_idx)',Data{3}(:,param_idx)'];
GroupNames = [string(repmat({groupname},1,size(Data{1},1))),...
                 string(repmat({groupname},1,size(Data{2},1))),...
                 string(repmat({groupname},1,size(Data{3},1)))];
ConditionNames = [string(repmat({'NoChange'},1,size(Data{1},1))),...
                 string(repmat({'NoDistalCue'},1,size(Data{2},1))),...
                 string(repmat({'NoOpticFlow'},1,size(Data{3},1)))];    
%remove NaN rows. There exsits nan rows because trialNums<paramNums of specific participant and specific condition
nonan_idx = ~isnan(Y);
Y = Y(nonan_idx);
GroupNames = GroupNames(nonan_idx);
ConditionNames = ConditionNames(nonan_idx);

end