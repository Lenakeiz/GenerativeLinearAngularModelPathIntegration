function TTest_NominalValue(groupData, nominalValue, groupName, tailtype)
%% TTestNominalValue
% Andrea Castegnaro, UCL, 2023 uceeaca@ucl.ac.uk
% Calculate ttest against a nominal value. Include possibility to decide
% which direction. Print the results.
% ===================================================================================

% Perform one-sample t-test
[h,p,ci,stats] = ttest(groupData, nominalValue, 'Tail', tailtype);

% Extract t-value and degrees of freedom
t_value = stats.tstat;
df = stats.df;

% Display the results
disp(['Results for group: ', groupName])
disp(['t(' num2str(df) ') = ' num2str(t_value) ', p = ' num2str(p)])
end