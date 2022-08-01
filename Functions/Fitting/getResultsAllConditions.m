%% A function for getting Results from All Conditions
function Results = getResultsAllConditions(TransformedData, config)

    if(config.useTrialFilter)
        %get the estimated parameters, X, DX, Theta, IC for all trial conditions for each group of data.
        for TRIAL_FILTER=1:3
            disp(['%%%%%%%%%%%%%%% PERFORMING FITTING: CONDITION ', num2str(TRIAL_FILTER), ' %%%%%%%%%%%%%%%']);
            config.TrialFilter                    =       TRIAL_FILTER;
            ModelFitResults                       =       PerformGroupFit(TransformedData, config);
            Results.estimatedParams{TRIAL_FILTER} =       ModelFitResults.estimatedParams; 
            Results.X{TRIAL_FILTER}               =       ModelFitResults.X;
            Results.DX{TRIAL_FILTER}              =       ModelFitResults.DX;      
            Results.THETADX{TRIAL_FILTER}         =       ModelFitResults.THETADX;
            Results.IC{TRIAL_FILTER}              =       ModelFitResults.IC;
            Results.DistErr{TRIAL_FILTER}         =       ModelFitResults.DistErr;
            Results.AngleErr{TRIAL_FILTER}        =       ModelFitResults.AngleErr;
            Results.PropDistErr{TRIAL_FILTER}     =       ModelFitResults.PropDistErr;
            Results.PropAngErr{TRIAL_FILTER}      =       ModelFitResults.PropAngErr;
            Results.LocationErr{TRIAL_FILTER}     =       ModelFitResults.LocationErr;
            Results.flagOoB{TRIAL_FILTER}         =       ModelFitResults.flagOoB;
            Results.L1Dur{TRIAL_FILTER}           =       ModelFitResults.L1Dur;
            Results.L2Dur{TRIAL_FILTER}           =       ModelFitResults.L2Dur;
        end
    else
        config.TrialFilter                        =       0;
        disp('%%%%%%%%%%%%%%% PERFORMING FITTING: MERGED CONDITIONS %%%%%%%%%%%%%%%');
        Results = PerformGroupFit(TransformedData, config);
    end
end