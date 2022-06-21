%% A function for getting Results from All Conditions
function [AllParams, AllX, AllDX, AllTheta, AllDistErr, AllAngleErr, AllFlagOoB, AllIC] = getResultsAllConditions(TransformedData, config)
    %get the estimated parameters, X, DX, Theta, IC for all trial conditions for each group of data.
    for TRIAL_FILTER=1:3
        disp('%%%%%%%%%%%%%%% PERFORMING FITTING %%%%%%%%%%%%%%%');
        config.TrialFilter              =       TRIAL_FILTER;
        Results                         =       PerformGroupFit(TransformedData, config);
        AllParams{TRIAL_FILTER}         =       Results.estimatedParams; 
        AllX{TRIAL_FILTER}              =       Results.X;
        AllDX{TRIAL_FILTER}             =       Results.DX;      
        AllTheta{TRIAL_FILTER}          =       Results.THETADX;
        AllIC{TRIAL_FILTER}             =       Results.IC;
        AllDistErr{TRIAL_FILTER}        =       Results.DistErr;
        AllAngleErr{TRIAL_FILTER}       =       Results.AngleErr;
        AllFlagOoB{TRIAL_FILTER}        =       Results.flagOoB;
    end
end    