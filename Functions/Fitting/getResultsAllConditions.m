%% A function for getting Results from All Conditions
function Results = getResultsAllConditions(TransformedData, config)
    %get the estimated parameters, X, DX, Theta, IC for all trial conditions for each group of data.
    for TRIAL_FILTER=1:3
        disp('%%%%%%%%%%%%%%% PERFORMING FITTING %%%%%%%%%%%%%%%');
        config.TrialFilter                    =       TRIAL_FILTER;
        ModelFitResults                       =       PerformGroupFit(TransformedData, config);
        Results{TRIAL_FILTER}.estimatedParams =       ModelFitResults.estimatedParams; 
        Results{TRIAL_FILTER}.X               =       ModelFitResults.X;
        Results{TRIAL_FILTER}.DX              =       ModelFitResults.DX;      
        Results{TRIAL_FILTER}.THETADX         =       ModelFitResults.THETADX;
        Results{TRIAL_FILTER}.IC              =       ModelFitResults.IC;
        Results{TRIAL_FILTER}.DistErr         =       ModelFitResults.DistErr;
        Results{TRIAL_FILTER}.AngleErr        =       ModelFitResults.AngleErr;
        Results{TRIAL_FILTER}.PropDistErr     =       ModelFitResults.PropDistErr;
        Results{TRIAL_FILTER}.PropAngErr      =       ModelFitResults.PropAngErr;
        Results{TRIAL_FILTER}.flagOoB         =       ModelFitResults.flagOoB;
    end
end