%% A function for getting Results from All Conditions
function [AllParams, AllX, AllDX, AllTheta, AllICDist] = getResultsAllConditions(TransformedData, config)
    %get the estimated parameters, X, DX, Theta, IC for all trial conditions for each group of data.
    AllParams = cell(0); AllX = cell(0);AllDX = cell(0); AllTheta = cell(0); AllICDist = cell(0);
    for TRIAL_FILTER=1:3
        config.TrialFilter = TRIAL_FILTER;
        tic
        disp('%%%%%%%%%%%%%%% PERFORMING FITTING %%%%%%%%%%%%%%%');
        Results = PerformGroupFit(TransformedData, config);
        AllParams{TRIAL_FILTER} = Results.estimatedParams; 
        AllX{TRIAL_FILTER} = Results.X;
        AllDX{TRIAL_FILTER}=Results.DX;      
        AllTheta{TRIAL_FILTER}=Results.THETADX;
        AllICDist{TRIAL_FILTER}=Results.ICDist;
        AllICAng{TRIAL_FILTER}=Results.ICDist;
        toc
    end
end    