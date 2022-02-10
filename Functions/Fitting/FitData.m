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
%% set model configurations

%set lower bound and up bound
%      1-gamma    2-G3    3-g2     4-g3     5-b      6-sigma    7-nu
lb  = [0.5,       0.1,    0.5,     0,       0,         0.1,       0.1];
ub  = [1.5,       1.0,    2.0,     1.0,     2*pi,      2.0,       100.0]; 

%set equality constriants
Aeq = zeros(7,7); beq=zeros(1,7);
%set the equality constriant of G3=1, this is common across models. G3
%controls the sacling of the mental triangle
Aeq(1,2)=1; beq(1)=1;

if Model_Name=="base"
    %do nothing
elseif Model_Name=="base_add_G3"
    %remove constraints on G3
    Aeq(1,2)=0; beq(1)=0;   
elseif Model_Name=="set_k3_0"
    %set equalities constriant of k3=0
    Aeq(2,5)=1;beq(2)=0;   
elseif Model_Name=="set_g3_1_k3_0"
    %set equalities constriant of g3=1 k3=0
    Aeq(2,4)=1;beq(2)=1;
    Aeq(3,5)=1;beq(3)=0;
elseif Model_Name=="set_g2_1"
    %set equalities constriant of g2=1
    Aeq(2,3)=1;beq(2)=1;          
elseif Model_Name=="set_same_sg"
    %set equalities constriant of g2=g3
    Aeq(2,3)=1;Aeq(2,4)=-1;beq(2)=0;  
    %set equalities constriant of k3=0
    Aeq(3,5)=1;beq(3)=0;  
elseif Model_Name=="allo" | Model_Name=="allo_weber"
    Aeq = eye(7,7); beq=ones(1,7); %identity matrix for Aeq
    Aeq(5,5)=1; beq(5)=0; %b=0,e.g., no bias
    Aeq(6,6)=0; beq(6)=0; %free sigma
else
    error("Please set the correct name of model!");
end

%calculate the likelihood function
if Model_Name=="allocentric"
    useweber=false;
    estFnc = @(FP) Estimate_allocentric(FP(6),DX,THETAX,X,useweber);
elseif Model_Name=="allocentric_weber"
    useweber=true;
    estFnc = @(FP) Estimate_allocentric(FP(6),DX,THETAX,X,useweber);   
else
    estFnc = @(FP) Estimate(FP(1),FP(2),FP(3),FP(4),FP(5),FP(6),FP(7),DX,THETAX,X,flagOoB);
end

%% parameter fitting

%Generating random start in the range
FitParams0 = (ub - lb)'.*rand(size(lb,2),1) + lb';
disp("Initial parameters: "+num2str(FitParams0'));
%setting optimization options
optim_options = optimoptions(@fmincon, 'Algorithm','sqp', 'MaxFunctionEvaluations',1e4);
%setting nonlinear constriants
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

