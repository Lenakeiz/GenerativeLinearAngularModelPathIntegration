function [Y, GroupNames, ConditionNames, GenderNames]=GroupAndRemoveNaN_3way_OnData(Data, groupname, gender)
%Group all the data from all groups and all conditions into a long
%numeric vector for further three-way anova analysis 
%and also output the names 

Y = [Data(:,1)',Data(:,2)',Data(:,3)'];
GroupNames = [string(repmat({groupname},1,size(Data,1))),...
                 string(repmat({groupname},1,size(Data,1))),...
                 string(repmat({groupname},1,size(Data,1)))];
ConditionNames = [string(repmat({'NoChange'},1,size(Data,1))),...
                 string(repmat({'NoDistalCue'},1,size(Data,1))),...
                 string(repmat({'NoOpticFlow'},1,size(Data,1)))];  
GenderNames = [gender', gender', gender'];

%remove NaN rows. There exsits nan rows because trialNums<paramNums of specific participant and specific condition
nonan_idx       =       ~isnan(Y);
Y               =       Y(nonan_idx);
GroupNames      =       GroupNames(nonan_idx);
ConditionNames  =       ConditionNames(nonan_idx);
GenderNames     =       GenderNames(nonan_idx);

end