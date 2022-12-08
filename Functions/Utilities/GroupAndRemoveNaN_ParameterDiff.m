function [Y, GroupNames, GenderNames]=GroupAndRemoveNaN_ParameterDiff(Data,param_idx, groupname,gender)
%Group all the data from all groups and all conditions into a long
%numeric vector for further two-way anova analysis 
%and also output the names 

Y = Data(:,param_idx)';
GroupNames = string(repmat({groupname},1,size(Data,1)));
GenderNames = gender';
%remove NaN rows. There exsits nan rows because trialNums<paramNums of specific participant and specific condition
nonan_idx = ~isnan(Y);
Y = Y(nonan_idx);
GroupNames = GroupNames(nonan_idx);
GenderNames= GenderNames(nonan_idx);

end