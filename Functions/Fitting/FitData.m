function [FitAllParams, ICDist, ICAng] = FitData(Input,config)
%FITDATA Function to fit the data on a single participant

%load configurations necessary for the script
Model_Name      =   config.ModelName;
useglobalsearch =   config.UseGlobalSearch;

if Model_Name=="AlloModel" 
    %set parameter lower bound and up bound
    %      1-sigma 
    lb  = [0.1];
    ub  = [2.0]; 
    estFnc = @(FP) EstimateAllo(FP,Input, config); 

%% ConstSpeedModelDistAngleGain model
elseif Model_Name=="ConstSpeedModelDistAngleGain"

    % distance parameter fitting

    %set parameter lower bound and up bound
    %     1, beta    2, sigma      
    lb_dist  = [-1.0,     0.1];
    ub_dist  = [1.0,      2.0]; 

    %calculate the likelihood function of distance error
    estFncDist = @(FP) EstimateConstSpeed_Dist(FP(1), FP(2), Input, config); 

    %Generating random start in the range
    FitDistParams0 = (ub_dist - lb_dist)'.*rand(size(lb_dist,2),1) + lb_dist';
    disp("Initial distance parameters: "+num2str(FitDistParams0'));
    %setting optimization options
    optim_options = optimoptions(@fmincon, 'Algorithm','sqp', 'MaxFunctionEvaluations',1e4);
    
    if useglobalsearch == true
        %finding the best local minima with globalsearch
        problem = createOptimProblem('fmincon', ...
                                     'objective', estFncDist, ...
                                     'x0', FitDistParams0, ...
                                     'Aineq',[],'bineq',[], ...
                                     'Aeq',[],'beq',[], ...
                                     'lb',lb_dist,'ub',ub_dist, ...
                                     'nonlcon', [], ...
                                     'options',optim_options);
        gs = GlobalSearch('NumTrialPoints', 1000);
        [FitDistParams,negloglikelihood_Dist] = run(gs,problem);
    else
        %find a local minima with fmincon
        [FitDistParams, negloglikelihood_Dist] = fmincon(estFncDist,FitDistParams0,[],[],[],[],lb_dist,ub_dist,nonlcon,optim_options);
    end
    
    disp("Fitted distance parameters: "+num2str(FitDistParams'));
    disp("Distance Negative LogLikelihood=" + num2str(negloglikelihood_Dist));
    disp(" ");disp(" ");disp(" ");

    % Calculate the Bayesian Inference Criterion for each model
    sampleSize = size(Input.X,2); %the number of observations, i.e., the sample size
    loglikelihood = -negloglikelihood_Dist; %the loglikelihood derived from fitting different models 
    numParams = length(lb_dist);
    [aic, bic] = aicbic(loglikelihood, numParams, sampleSize, 'Normalize',false);
    ICDist.aic = aic;
    ICDist.bic = bic;
    ICDist.negll = negloglikelihood_Dist;
    ICDist.likelihood = exp(loglikelihood);    

    % angle parameter fitting
    
    beta = FitDistParams(1); %from distance estimation

    %           1-g3        2-nu
    lb_angle  = [0,         0.1];
    ub_angle  = [2.0,       100.0]; 

    %calculate the likelihood function of angle error
    estFncAng = @(FP) EstimateConstSpeed_AngleGain(FP(1), FP(2), beta, Input, config); 

    %Generating random start in the range
    FitAngParams0 = (ub_angle - lb_angle)'.*rand(size(lb_angle,2),1) + lb_angle';
    disp("Initial angle parameters: "+num2str(FitAngParams0'));
    %setting optimization options
    optim_options = optimoptions(@fmincon, 'Algorithm','sqp', 'MaxFunctionEvaluations',1e4);
    
    if useglobalsearch == true
        %finding the best local minima with globalsearch
        problem = createOptimProblem('fmincon', ...
                                     'objective', estFncAng, ...
                                     'x0', FitAngParams0, ...
                                     'Aineq',[],'bineq',[], ...
                                     'Aeq',[],'beq',[], ...
                                     'lb',lb_angle,'ub',ub_angle, ...
                                     'nonlcon', [], ...
                                     'options',optim_options);
        gs = GlobalSearch('NumTrialPoints', 1000);
        [FitAngParams,negloglikelihood_Ang] = run(gs,problem);
    else
        %find a local minima with fmincon
        [FitAngParams, negloglikelihood_Ang] = fmincon(estFncAng,FitAngParams0,[],[],[],[],lb_angle,ub_angle,nonlcon,optim_options);
    end
    
    disp("Fitted angle parameters: "+num2str(FitAngParams'));
    disp("Angle Negative LogLikelihood=" + num2str(negloglikelihood_Ang));
    disp(" ");disp(" ");disp(" ");    

    % Calculate the Bayesian Inference Criterion for each model
    sampleSize = size(Input.X,2); %the number of observations, i.e., the sample size
    loglikelihood = -negloglikelihood_Ang; %the loglikelihood derived from fitting different models 
    numParams = length(lb_angle);
    [aic, bic] = aicbic(loglikelihood, numParams, sampleSize, 'Normalize',false);
    ICAng.aic = aic;
    ICAng.bic = bic;
    ICAng.negll = negloglikelihood_Ang;
    ICAng.likelihood = exp(loglikelihood); 

    %
    FitAllParams = [FitDistParams; FitAngParams];

%% ConstSpeedModelDistAngleRGb model
elseif Model_Name=="ConstSpeedModelDistAngleRGb"

    %% distance parameter fitting

    %set parameter lower bound and up bound
    %     1, beta    2, sigma      
    lb_dist  = [-1.0,     0.1];
    ub_dist  = [1.0,      2.0]; 

    %calculate the likelihood function of distance error
    estFncDist = @(FP) EstimateConstSpeed_Dist(FP(1), FP(2), Input, config); 

    %Generating random start in the range
    FitDistParams0 = (ub_dist - lb_dist)'.*rand(size(lb_dist,2),1) + lb_dist';
    disp("Initial distance parameters: "+num2str(FitDistParams0'));
    %setting optimization options
    optim_options = optimoptions(@fmincon, 'Algorithm','sqp', 'MaxFunctionEvaluations',1e4);
    
    if useglobalsearch == true
        %finding the best local minima with globalsearch
        problem = createOptimProblem('fmincon', ...
                                     'objective', estFncDist, ...
                                     'x0', FitDistParams0, ...
                                     'Aineq',[],'bineq',[], ...
                                     'Aeq',[],'beq',[], ...
                                     'lb',lb_dist,'ub',ub_dist, ...
                                     'nonlcon', [], ...
                                     'options',optim_options);
        gs = GlobalSearch('NumTrialPoints', 1000);
        [FitDistParams,negloglikelihood_Dist] = run(gs,problem);
    else
        %find a local minima with fmincon
        [FitDistParams, negloglikelihood_Dist] = fmincon(estFncDist,FitDistParams0,[],[],[],[],lb_dist,ub_dist,nonlcon,optim_options);
    end
    
    disp("Fitted distance parameters: "+num2str(FitDistParams'));
    disp("Distance Negative LogLikelihood=" + num2str(negloglikelihood_Dist));
    disp(" ");disp(" ");disp(" ");

    % Calculate the Bayesian Inference Criterion for each model
    sampleSize = size(Input.X,2); %the number of observations, i.e., the sample size
    loglikelihood = -negloglikelihood_Dist; %the loglikelihood derived from fitting different models 
    numParams = length(lb_dist);
    [aic, bic] = aicbic(loglikelihood, numParams, sampleSize, 'Normalize',false);
    ICDist.aic = aic;
    ICDist.bic = bic;
    ICDist.negll = negloglikelihood_Dist;
    ICDist.likelihood = exp(loglikelihood);    

    %% angle parameter fitting
    
    beta = FitDistParams(1); %from distance estimation

    %           1-g3    2-b     3-nu
    lb_angle  = [0,     0,      0.1];
    ub_angle  = [2.0,   2*pi,   100.0]; 

    %calculate the likelihood function of angle error
    estFncAng = @(FP) EstimateConstSpeed_AngleRGb(FP(1), FP(2), FP(3), beta, Input, config); 

    %Generating random start in the range
    FitAngParams0 = (ub_angle - lb_angle)'.*rand(size(lb_angle,2),1) + lb_angle';
    disp("Initial angle parameters: "+num2str(FitAngParams0'));
    %setting optimization options
    optim_options = optimoptions(@fmincon, 'Algorithm','sqp', 'MaxFunctionEvaluations',1e4);
    
    if useglobalsearch == true
        %finding the best local minima with globalsearch
        problem = createOptimProblem('fmincon', ...
                                     'objective', estFncAng, ...
                                     'x0', FitAngParams0, ...
                                     'Aineq',[],'bineq',[], ...
                                     'Aeq',[],'beq',[], ...
                                     'lb',lb_angle,'ub',ub_angle, ...
                                     'nonlcon', [], ...
                                     'options',optim_options);
        gs = GlobalSearch('NumTrialPoints', 1000);
        [FitAngParams,negloglikelihood_Ang] = run(gs,problem);
    else
        %find a local minima with fmincon
        [FitAngParams, negloglikelihood_Ang] = fmincon(estFncAng,FitAngParams0,[],[],[],[],lb_angle,ub_angle,nonlcon,optim_options);
    end
    
    disp("Fitted angle parameters: "+num2str(FitAngParams'));
    disp("Angle Negative LogLikelihood=" + num2str(negloglikelihood_Ang));
    disp(" ");disp(" ");disp(" ");    

    % Calculate the Bayesian Inference Criterion for each model
    sampleSize = size(Input.X,2); %the number of observations, i.e., the sample size
    loglikelihood = -negloglikelihood_Ang; %the loglikelihood derived from fitting different models 
    numParams = length(lb_angle);
    [aic, bic] = aicbic(loglikelihood, numParams, sampleSize, 'Normalize',false);
    ICAng.aic = aic;
    ICAng.bic = bic;
    ICAng.negll = negloglikelihood_Ang;
    ICAng.likelihood = exp(loglikelihood); 

    %
    FitAllParams = [FitDistParams; FitAngParams];

elseif Model_Name=="ConstSpeedModel5Params"
    %set parameter lower bound and up bound
    %     1, beta        2-g3     3-b      4-sigma      5-nu
    lb  = [-1.0,         0,       0,      0.1,         0.1];
    ub  = [1.0,          2.0,    2*pi,    2.0,         100.0];    

    %set equality constriants
    Aeq         =   [];         beq     =   [];

    %calculate the likelihood function
    estFnc = @(FP) EstimateConstSpeed5Params(FP(1),FP(2),FP(3),FP(4),FP(5), Input, config); 

elseif Model_Name == "ConstSpeedModelwith_g2"
    %set parameter lower bound and up bound
    %     1, beta        2-g2     3-g3      4-sigma      5-nu
    lb  = [-1.0,         0,       0,        0.1,         0.1];
    ub  = [1.0,          3,       3,        2.0,         3.14]; 

    %calculate the likelihood function
    estFnc = @(FP) EstimateConstSpeedwith_g2(FP(1),FP(2),FP(3),FP(4),FP(5), Input, config);  

    %Generating random start in the range
    FitParams0 = (ub - lb)'.*rand(size(lb,2),1) + lb';
    disp("Initial parameters: "+num2str(FitParams0'));
    %setting optimization options
    optim_options = optimoptions(@fmincon, 'Algorithm','sqp', 'MaxFunctionEvaluations',1e4);

    %finding the best local minima with globalsearch
    problem = createOptimProblem('fmincon','objective', estFnc,'x0',FitParams0, ...
                                 'Aineq',[],'bineq',[], ...
                                 'Aeq',[],'beq',[], ...
                                 'lb',lb,'ub',ub, ...
                                 'nonlcon', [], ...
                                 'options',optim_options);
    gs = GlobalSearch('NumTrialPoints', 1000);
    [FitAllParams,negloglikelihood] = run(gs,problem);

    disp("Fitted parameters: "+num2str(FitAllParams'));
    disp("Negative LogLikelihood=" + num2str(negloglikelihood));
    disp(" ");disp(" ");disp(" ");
    ICDist = 0;
    ICAng  = 0;

elseif Model_Name == "ConstSpeedModelwith_g2_k3"
    %set parameter lower bound and up bound
    %     1, beta        2-g2     3-g3    4-k3    5-sigma      6-nu
    lb  = [-1.0,         0,       0,      0,      0.0,         0.0];
    ub  = [1.0,          3,       3,      pi,     4.0,         pi]; 

    %calculate the likelihood function
    estFnc = @(FP) EstimateConstSpeedwith_g2_k3(FP(1),FP(2),FP(3),FP(4),FP(5), FP(6), Input, config);  

    %Generating random start in the range
    FitParams0 = (ub - lb)'.*rand(size(lb,2),1) + lb';
    disp("Initial parameters: "+num2str(FitParams0'));
    %setting optimization options
    optim_options = optimoptions(@fmincon, 'Algorithm','sqp', 'MaxFunctionEvaluations',1e4);

    %finding the best local minima with globalsearch
    problem = createOptimProblem('fmincon','objective', estFnc,'x0',FitParams0, ...
                                 'Aineq',[],'bineq',[], ...
                                 'Aeq',[],'beq',[], ...
                                 'lb',lb,'ub',ub, ...
                                 'nonlcon', [], ...
                                 'options',optim_options);
    gs = GlobalSearch('NumTrialPoints', 1000);
    [FitAllParams,negloglikelihood] = run(gs,problem);

    disp("Fitted parameters: "+num2str(FitAllParams'));
    disp("Negative LogLikelihood=" + num2str(negloglikelihood));
    disp(" ");disp(" ");disp(" ");
    ICDist = 0;
    ICAng  = 0;

elseif Model_Name == "ConstSpeedModelwith_g2_RGmean"
    %set parameter lower bound and up bound
    %     1, beta        2-g2     3-g3    4-sigma      5-nu
    lb  = [-1.0,         0,       0,      0.0,         0.0];
    ub  = [1.0,          3,       3,      4.0,         pi]; 

    %calculate the likelihood function
    estFnc = @(FP) EstimateConstSpeedwith_g2_RGmean(FP(1),FP(2),FP(3),FP(4),FP(5), Input, config);  

    %Generating random start in the range
    FitParams0 = (ub - lb)'.*rand(size(lb,2),1) + lb';
    disp("Initial parameters: "+num2str(FitParams0'));
    %setting optimization options
    optim_options = optimoptions(@fmincon, 'Algorithm','sqp', 'MaxFunctionEvaluations',1e4);

    %finding the best local minima with globalsearch
    problem = createOptimProblem('fmincon','objective', estFnc,'x0',FitParams0, ...
                                 'Aineq',[],'bineq',[], ...
                                 'Aeq',[],'beq',[], ...
                                 'lb',lb,'ub',ub, ...
                                 'nonlcon', [], ...
                                 'options',optim_options);
    gs = GlobalSearch('NumTrialPoints', 1000);
    [FitAllParams,negloglikelihood] = run(gs,problem);

    disp("Fitted parameters: "+num2str(FitAllParams'));
    disp("Negative LogLikelihood=" + num2str(negloglikelihood));
    disp(" ");disp(" ");disp(" ");
    ICDist = 0;
    ICAng  = 0;

elseif Model_Name == "IntSpeedModel"
    %set parameter lower bound and up bound
    %     1, beta    2-G3     3-g2     4-g3     5-b      6-sigma      7-nu
    lb  = [-1.0,     0.5,     0.5,     0.2,       0,       0.1,         0.1];
    ub  = [1.0,      2.0,     2.0,     2.0,     2*pi,     2.0,         100.0];    

    %set equality constriants
    Aeq         =   zeros(7,7);         beq     =   zeros(1,7);
    Aeq(2,2)    =   1;                  beq(2)  =   1;             %G3=1
    Aeq(3,3)    =   1;                  beq(3)  =   1;             %g2=1

    if config.subtype == "DistAng_RGb" 
        %do nothing
    elseif config.subtype == "DistAng_RGmean"
        Aeq(5,5) = 1;  beq(5) = 0; %b=0
    else
        error("Please set the correct subtype!");
    end

    %calculate the likelihood function
    estFnc = @(FP) EstimateIntSpeed(FP(1),FP(2),FP(3),FP(4),FP(5),FP(6),FP(7),Input, config);

elseif Model_Name=="GammaModel" 
    %set parameter lower bound and up bound
    %      1-gamma    2-G3    3-g2     4-g3     5-b      6-sigma    7-nu
    lb  = [0.5,       0.1,    0.5,     0,       0,         0.1,       0.1];
    ub  = [1.5,       1.0,    2.0,     2.0,     pi,      2.0,       100.0]; 

    %set equality constriants
    Aeq         =   zeros(7,7);         beq     =   zeros(1,7);
    Aeq(2,2)    =   1;                  beq(2)  =   1;             %G3=1
    Aeq(3,3)    =   1;                  beq(3)  =   1;             %g2=1

    if config.subtype == "DistAng_RGb" 
        %do nothing
    elseif config.subtype == "DistAng_RGmean"
        Aeq(5,5) = 1;  beq(5) = 0; %b=0
    else
        error("Please set the correct subtype!");
    end

    %calculate the likelihood function
    estFnc = @(FP) EstimateGamma(FP(1),FP(2),FP(3),FP(4),FP(5),FP(6),FP(7), Input, config);

elseif Model_Name == "G1G2Model"
    %set model configurations
    %set lower bound and up bound
    %      1-G1     2-G2    3-G3    4-g2   5-g3   6-b    7-sigma    8-nu
    lb  = [0,       0,      0.1,    0.5,   0,     0,         0.1,       0.1];
    ub  = [1.5,     1.5,    1.0,    2.0,   2.0,   pi,        2.0,       100.0];

    %set equality constriants
    Aeq         =   zeros(8,8);         beq     =   zeros(1,8);
    Aeq(3,3)    =   1;                  beq(3)  =   1;             %G3=1
    Aeq(4,4)    =   1;                  beq(4)  =   1;             %g2=1

    if config.subtype == "DistAng_RGb" 
        %do nothing
    elseif config.subtype == "DistAng_RGmean"
        Aeq(6,6) = 1;  beq(6) = 0; %b=0
    else
        error("Please set the correct subtype!");
    end

    %calculate the likelihood function
    estFnc = @(FP) EstimateG1G2(FP(1),FP(2),FP(3),FP(4),FP(5),FP(6),FP(7),FP(8), Input, config);

else
    error("Please set the correct name of model!");
end

end

