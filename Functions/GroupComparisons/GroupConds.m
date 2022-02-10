function [Y, GroupNames, ConditionNames]=GroupConds(Data,groupname)
%Group all the data from five groups and three conditions into a long
%numeric vector for further two-way anova analysis 
%and also output the factor names 


Y = [Data{1},Data{2},Data{3}];
GroupNames = [string(repmat({groupname},1,size(Data{1},2))),...
                 string(repmat({groupname},1,size(Data{2},2))),...
                 string(repmat({groupname},1,size(Data{3},2)))];
ConditionNames = [string(repmat({'NoChange'},1,size(Data{1},2))),...
                 string(repmat({'NoDistalCue'},1,size(Data{2},2))),...
                 string(repmat({'NoOpticFlow'},1,size(Data{3},2)))];    

end