function [FitParams, IC] = FitData(DX,THETAX,X,OoBLen,flagOoB,config)
%FITDATA Function to fit the data on a single participant
%   DX is a cell structure containing the segment of each trial
%   THETAX is the turning angle (wrong at the moment)
%   X is the data points
%   OoBLen is the Out-of-Boundary length with 0 for non-oob trials, >0 for oob trials.
%   flagOoB is the flag of Out-of-Boundary trials, with 0 non-OoB trials 1 those OoB trials
%   the noise for each trial

%load configurations necessary for the script
Model_Name = config.ModelName;
numParams = config.NumParams;
useglobalsearch = config.UseGlobalSearch;
ifallo=false; %only true when the model is a simple allocentric generative model
useweber=false;%only true when use weber law in simple generative models

if Model_Name == "G1G2"
    %set model configurations
    %set lower bound and up bound
    %      1-G1     2-G2    3-G3    4-g2   5-g3   6-b    7-sigma    8-nu
    lb  = [0.1,     0.1,    0.1,    0.5,   0,     0,         0.1,       0.1];
    ub  = [1.0,     1.0,    1.0,    2.0,   1.0,   2*pi,      2.0,     100.0];

    %set equality constriants
    Aeq = zeros(8,8); beq=zeros(1,8); 
    
    if Model_Name=="G1G2"
        Aeq(3,3)=1; beq(3)=1;%G3=1
        Aeq(4,4)=1; beq(4)=1;%g2=1        
    else
        error("Please set the correct name of model!");
    end    
    %calculate the likelihood function
    estFnc = @(FP) EstimateG1G2(FP(1),FP(2),FP(3),FP(4),FP(5),FP(6),FP(7),FP(8),DX,THETAX,X);

else %gamma based model
    %set model configurations
    %set lower bound and up bound
    %      1-gamma    2-G3    3-g2     4-g3     5-b      6-sigma    7-nu
    lb  = [0.5,       0.1,    0.5,     0,       0,         0.1,       0.1];
    ub  = [1.5,       1.0,    2.0,     1.0,     2*pi,      2.0,       100.0]; 

    %set equality constriants
    Aeq = zeros(7,7); beq=zeros(1,7);
    
    if Model_Name=="BaseModel"
        Aeq(2,2)=1; beq(2)=1;%G3=1
        Aeq(3,3)=1; beq(3)=1;%g2=1
    elseif Model_Name=="DistErrModel"
        Aeq(2,2)=1; beq(2)=1;%G3=1
        Aeq(3,3)=1; beq(3)=1;%g2=1 
        Aeq(4,4)=1; beq(4)=1;%g3=1 
        Aeq(5,5)=1; beq(5)=0;%b=0 
    elseif Model_Name=="AngleErrModel"
        Aeq(1,1)=1; beq(1)=1;%gamma=1
        Aeq(2,2)=1; beq(2)=1;%G3=1
        Aeq(3,3)=1; beq(3)=1;%g2=1             
    elseif Model_Name=="Allo" | Model_Name=="AlloWeber"
        Aeq(1,1)=1; beq(1)=1;%gamma=1
        Aeq(2,2)=1; beq(2)=1;%G3=1
        Aeq(3,3)=1; beq(3)=1;%g2=1   
        Aeq(4,4)=1; beq(4)=1;%g3=1 
        Aeq(5,5)=1; beq(5)=0;%b=0 
        Aeq(7,7)=1; beq(7)=0;%b=0 
        ifallo=true;
        if Model_Name=="AlloWeber"
            useweber=true;
        end
    elseif Model_Name=="Ego" | Model_Name=="EgoWeber"
        Aeq(1,1)=1; beq(1)=1;%gamma=1
        Aeq(2,2)=1; beq(2)=1;%G3=1
        Aeq(3,3)=1; beq(3)=1;%g2=1   
        Aeq(4,4)=1; beq(4)=1;%g3=1 
        Aeq(5,5)=1; beq(5)=0;%b=0 
        if Model_Name=="EgoWeber"
            useweber=true;
        end    
    else
        error("Please set the correct name of model!");
    end
    
    %calculate the likelihood function
    estFnc = @(FP) Estimate(FP(1),FP(2),FP(3),FP(4),FP(5),FP(6),FP(7),DX,THETAX,X,ifallo,useweber);
end

%% parameter fitting
%Generating random start in the range
FitParams0 = (ub - lb)'.*rand(size(lb,2),1) + lb';
disp("Initial parameters: "+num2str(FitParams0'));
%setting optimization options
optim_options = optimoptions(@fmincon, 'Algorithm','sqp', 'MaxFunctionEvaluations',1e4);

%setting nonlinear constriants We don't have to constraint here coz
%Theta3Prime will always in the range of (0,2pi)
%nonlcon = @(FP) Theta3PrimeCon(FP(1),FP(2),FP(3),FP(4),FP(5),DX,THETAX,X);
nonlcon = [];

if useglobalsearch == true
    %finding the best local minima with globalsearch
    problem = createOptimProblem('fmincon','objective', estFnc,'x0',FitParams0, ...
                                 'Aineq',[],'bineq',[], ...
                                 'Aeq',Aeq,'beq',beq, ...
                                 'lb',lb,'ub',ub, ...
                                 'nonlcon', nonlcon, ...
                                 'options',optim_options);
    gs = GlobalSearch();
    [FitParams,negloglikelihood] = run(gs,problem);
else
    %find a local minima with fmincon
    [FitParams, negloglikelihood] = fmincon(estFnc,FitParams0,[],[],Aeq,beq,lb,ub,nonlcon,optim_options);
end

disp("Fitted parameters: "+num2str(FitParams'));
disp("Negative LogLikelihood=" + num2str(negloglikelihood));
disp(" ");
disp(" ");
disp(" ");

%% Calculate the Bayesian Inference Criterion for each model
sampleSize = size(X,2); %the number of observations, i.e., the sample size
loglikelihood = -negloglikelihood; %the loglikelihood derived from fitting different models 
[aic, bic] = aicbic(loglikelihood, numParams, sampleSize, 'Normalize',false);
IC.aic = aic;
IC.bic = bic;
IC.negll = negloglikelihood;
IC.likelihood = exp(loglikelihood);

end

