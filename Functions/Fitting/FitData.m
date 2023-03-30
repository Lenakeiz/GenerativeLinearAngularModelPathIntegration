function [FitAllParams, IC] = FitData(Input,config)
% Zilong Ji, UCL, zilong.ji@ucl.ac.uk
% Wrapper function 
% ===================================================================================

%load configurations necessary for the script
Model_Name      =   config.ModelName;

if Model_Name == "sigma_nu"      %simple egocentric noise model
    %set parameter lower bound and up bound
    %      1-sigma    2-nu
    lb  = [0.0,       0.0];
    ub  = [3.0,       pi]; 

    %defining likelihood function
    estFnc = @(FP) Estimate_sigma_nu(FP(1), FP(2), Input, config);

elseif Model_Name == "beta_k_sigma_nu"      % memory leaky decay + noise
    %set parameter lower bound and up bound
    %     1, beta     2-k    3-sigma    4-nu
    lb  = [-1.0,      0.5,     0.0,       0.0];
    ub  = [1.0,       2,       3.0,       pi]; 

    %defining likelihood function
    estFnc = @(FP) Estimate_beta_k_sigma_nu(FP(1), FP(2), FP(3), FP(4), Input, config);

elseif Model_Name == "g2_sigma_nu"      % memory leaky decay + noise
    %set parameter lower bound and up bound
    %     1-g2    2-sigma    3-nu
    lb  = [0.0,      0.0,       0.0];
    ub  = [3.0,      3.0,       pi]; 
    estFnc = @(FP) Estimate_g2_sigma_nu(FP(1), FP(2), FP(3), Input, config);

elseif Model_Name == "g3_sigma_nu"      %regressing to correct mean return angle
    %set parameter lower bound and up bound
    %      1-g3     2-sigma,     3-nu
    lb  = [0,       0.0,         0.0];
    ub  = [3,       3.0,         pi]; 

    %defining likelihood function
    estFnc = @(FP) Estimate_g3_sigma_nu(FP(1),FP(2),FP(3), Input, config);

elseif Model_Name == "m3_sigma_nu"      %regressing to correct mean return angle
    %set parameter lower bound and up bound
    %      1-m3     2-sigma,     3-nu
    lb  = [0,       0.0,         0.0];
    ub  = [3,       3.0,         pi]; 

    %defining likelihood function
    estFnc = @(FP) Estimate_m3_sigma_nu(FP(1),FP(2),FP(3), Input, config);

elseif Model_Name == "beta_k_g2_sigma_nu"      % memory leaky decay + noise
    %set parameter lower bound and up bound
    %     1, beta     2-k     3-g2    4-sigma    5-nu
    lb  = [-1.0,      0.5,    0.0,      0.0,       0.0];
    ub  = [1.0,       2.0,    3.0,      3.0,       pi]; 

    %defining likelihood function
    estFnc = @(FP) Estimate_beta_k_g2_sigma_nu(FP(1), FP(2), FP(3), FP(4), FP(5), Input, config);

elseif Model_Name == "beta_k_g3_sigma_nu"      % memory leaky decay + noise
    %set parameter lower bound and up bound
    %     1, beta     2-k     3-g3    4-sigma    5-nu
    lb  = [-1.0,      0.5,    0.0,      0.0,       0.0];
    ub  = [1.0,       2.0,    3.0,      3.0,       pi]; 

    %defining likelihood function
    estFnc = @(FP) Estimate_beta_k_g3_sigma_nu(FP(1), FP(2), FP(3), FP(4), FP(5), Input, config);

elseif Model_Name == "beta_k_m3_sigma_nu"      % memory leaky decay + noise
    %set parameter lower bound and up bound
    %     1, beta     2-k     3-m3    4-sigma    5-nu
    lb  = [-1.0,      0.5,    0.0,      0.0,       0.0];
    ub  = [1.0,       2.0,    3.0,      3.0,       pi]; 

    %defining likelihood function
    estFnc = @(FP) Estimate_beta_k_m3_sigma_nu(FP(1), FP(2), FP(3), FP(4), FP(5), Input, config);

elseif Model_Name == "g2_g3_sigma_nu"      %regressing to correct mean return angle
    %set parameter lower bound and up bound
    %      1-g2     2-g3     3-sigma,     4-nu
    lb  = [0,       0,       0.0,         0.0];
    ub  = [3,       3,       3.0,         pi]; 

    %defining likelihood function
    estFnc = @(FP) Estimate_g2_g3_sigma_nu(FP(1),FP(2),FP(3),FP(4), Input, config);

elseif Model_Name == "g2_m3_sigma_nu"      %regressing to correct mean return angle
    %set parameter lower bound and up bound
    %      1-g2     2-m3     3-sigma,     4-nu
    lb  = [0,       0,       0.0,         0.0];
    ub  = [3,       3,       3.0,         pi]; 

    %defining likelihood function
    estFnc = @(FP) Estimate_g2_m3_sigma_nu(FP(1),FP(2),FP(3),FP(4), Input, config);

elseif Model_Name == "g3_m3_sigma_nu"      %regressing to correct mean return angle
    %set parameter lower bound and up bound
    %      1-g3     2-m3     3-sigma,     4-nu
    lb  = [0,       0,       0.0,         0.0];
    ub  = [3,       3,       3.0,         pi]; 

    %defining likelihood function
    estFnc = @(FP) Estimate_g3_m3_sigma_nu(FP(1),FP(2),FP(3),FP(4), Input, config);    

elseif Model_Name == "beta_k_g2_g3_sigma_nu"      %regressing to correct mean return angle
    %set parameter lower bound and up bound
    %     1, beta       k         2-g2     3-g3     4-sigma      5-nu
    lb  = [-1.0,        0.5,         0,       0,       0.0,         0.0];
    ub  = [1.0,         2.0,         3,       3,       3.0,         pi];  
    %defining likelihood function
    estFnc = @(FP) Estimate_beta_k_g2_g3_sigma_nu(FP(1),FP(2),FP(3),FP(4),FP(5), FP(6),Input, config);  

elseif Model_Name == "beta_k_g2_m3_sigma_nu"      %regressing to correct mean return angle
    %set parameter lower bound and up bound
    %     1, beta       k         2-g2     3-m3     4-sigma      5-nu
    lb  = [-1.0,        0.5,      0,       0,       0.0,         0.0];
    ub  = [1.0,         2,       3,        3,       3.0,         pi];  
    %defining likelihood function
    estFnc = @(FP) Estimate_beta_k_g2_m3_sigma_nu(FP(1),FP(2),FP(3),FP(4),FP(5), FP(6),Input, config); 

elseif Model_Name == "beta_k_g3_m3_sigma_nu"      %regressing to correct mean return angle
    %set parameter lower bound and up bound
    %     1, beta       k         2-g3     3-m3     4-sigma      5-nu
    lb  = [-1.0,        0.5,      0,       0,       0.0,         0.0];
    ub  = [1.0,         2,        3,       3,       3.0,         pi];  
    %defining likelihood function
    estFnc = @(FP) Estimate_beta_k_g3_m3_sigma_nu(FP(1),FP(2),FP(3),FP(4),FP(5), FP(6),Input, config); 

elseif Model_Name == "g2_g3_m3_sigma_nu"      %regressing to correct mean return angle
    %set parameter lower bound and up bound
    %     1,g2        2-g3     3-m3     4-sigma      5-nu
    lb  = [0,         0,       0,       0.0,         0.0];
    ub  = [3,         3,       3,       3.0,         pi];  
    %defining likelihood function
    estFnc = @(FP) Estimate_g2_g3_m3_sigma_nu(FP(1),FP(2),FP(3),FP(4),FP(5), Input, config); 

elseif Model_Name == "beta_k_g2_g3_m3_sigma_nu"
    %set parameter lower bound and up bound
    %     1, beta     2-k      3-g2     4-g3    5-m3     6-sigma      7-nu
    lb  = [-1.0,      0.5,       0,       0,      0,       0.0,         0.0];
    ub  = [1.0,       2.0,       3,       3,      3,       3.0,         pi]; 

    %calculate the likelihood function
    estFnc = @(FP) Estimate_beta_k_g2_g3_m3_sigma_nu(FP(1),FP(2),FP(3),FP(4),FP(5),FP(6), FP(7), Input, config);    

elseif Model_Name == "beta1_beta2_g2_g3_sigma_nu"      %regressing to correct mean return angle
    %set parameter lower bound and up bound
    %     1, beta1       beta2         2-g2     3-g3     4-sigma      5-nu
    lb  = [0.2,           0.2,             0,       0,       0.0,         0.0];
    ub  = [2.5,           2.5,             3,       3,       3.0,         pi];  
    %defining likelihood function
    estFnc = @(FP) Estimate_beta1_beta2_g2_g3_sigma_nu(FP(1),FP(2),FP(3),FP(4),FP(5), FP(6),Input, config);  
else
    error("Please set the correct name of model!");
end

% Initializing parameters
FitParams0 = (ub - lb)'.*rand(size(lb,2),1) + lb';
disp("Initial parameters: "+num2str(FitParams0'));

% Setting optimization options
optim_options = optimoptions(@fmincon, 'MaxFunctionEvaluations',3e3);

% Finding the best local minima with globalsearch
problem = createOptimProblem('fmincon', ...
                             'objective', estFnc, ...
                             'x0',FitParams0, ...
                             'lb',lb,'ub',ub, ...
                             'options',optim_options);

gs = GlobalSearch('NumTrialPoints', 1000);
[FitAllParams,negloglikelihood] = run(gs,problem);

disp("Fitted parameters: "+num2str(FitAllParams'));
disp("Negative LogLikelihood=" + num2str(negloglikelihood));
disp(" ");disp(" ");disp(" ");

%% Calculate the Bayesian Inference Criterion for each model
sampleSize = size(Input.X,2); %the number of observations, i.e., the sample size

% Saving negative log likelihood
IC.negll = negloglikelihood;

% Saving likelihood
loglikelihood = -negloglikelihood; %the loglikelihood derived from fitting different models 
IC.likelihood = exp(loglikelihood);

% Saving AIC and BIC
[aic, bic, ~] = aicbic(loglikelihood, config.NumParams, sampleSize, 'Normalize',false);
IC.aic = aic;
IC.bic = bic;
end

