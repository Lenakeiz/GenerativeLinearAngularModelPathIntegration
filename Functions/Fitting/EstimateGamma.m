function [negloglikelihood] = EstimateGamma(gamma, G3, g2, g3, b, sigma, nu, Input, config)
%   EstimateGamma means a gamma discount factor on the walking distance,
%   with leg 1 discount more than leg 2
%   Args:
%       gamma is the discount factor
%       G3 is the gain of the length of the third leg 
%       g2 is the rotation angle from the direction of leg 2
%       g3 is the rotation angle from the direction of leg 3
%       b is the systematic bias in the execution error
%       sigma is the standard deviation for the Gaussian distribution of the return point
%       nu decribes the noise strength in the Von Mises distribution
%       Input contains all the data information for estimating, see PerformGroupFit for how it was generated
%       config: self-explained

%% information necessary for running parameter estimation
X        =   Input.X;
DX       =   Input.DX;
THETAX   =   Input.THETADX;
BadExecutions = Input.BadExecutionsTrial;
sampleSize = size(X,2);
negloglikelihood = 0;

%% find the correct mean return angle based on all trials 
Alphas = zeros(sampleSize,1);
for tr = 1:sampleSize
    %extract the physical data info
    l1      =   DX{tr}(1);
    l2      =   DX{tr}(2);
    theta2  =   THETAX{tr}(2);

    %calculate the correct return angle
    phy_p1  =   [l1,0];
    phy_p2  =   [l1+l2*cos(theta2),l2*sin(theta2)];
    vec1    =   phy_p2-phy_p1; vec2 = [0,0]-phy_p2;
    alpha   =   atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
    alpha   =   deg2rad(alpha);%transfer from degree to radians
    alpha   =   mod(alpha, 2*pi);  %wrap to (0,2pi)  
    Alphas(tr) = alpha;
end
mean_angle = mean(Alphas);

for tr = 1:sampleSize
    %% extract the physical data info
    l1      =   DX{tr}(1); 
    l2      =   DX{tr}(2); 
    l3      =   DX{tr}(3);
    theta2  =   THETAX{tr}(2); 
    theta3  =   THETAX{tr}(3); 

    %% whether to use weber's law to scaling the noise strength
    if config.useweber == true
        scale = sqrt(l1^2+l2^2); %noise scale
    else
        scale = 1; %noise scale
    end
    sigma_scaled = scale*sigma;
    nu_scaled = scale*nu;

    %mental point 1
    G1 = gamma^2*G3; 
    men_p1 = [G1*l1,0];

    theta2_prime = g2*theta2;

    %mental point 2
    G2 = gamma*G3;    
    men_p2 = [G1*l1+G2*l2*cos(theta2_prime),G2*l2*sin(theta2_prime)];

    %calculate length of mental vector 3
    h       =   norm(men_p2);
    %calculate turn angle of mental vector 3
    vec1    =   men_p2-men_p1; 
    vec2    =   [0,0]-men_p2;
    alpha   =   atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
    %transfer from degree to radians
    alpha   =   deg2rad(alpha);
    %wrap to (0,2pi)
    alpha   =   mod(alpha, 2*pi);

    %whether to regress to the mean correct return angle
    if contains("RGmean", config.subtype)
        theta3_prime = g3*alpha+(1-g3)*mean_angle;
    else 
        theta3_prime = g3*alpha+b;
    end    

    %also wrap theta3_prime to (0,2pi) (We don't have to coz of the 'cos' below, 
    %also coz of the nonconstriant (if there is) which restrict theta3_prime in [0,2pi]).
    theta3_prime = mod(theta3_prime, 2*pi);
    
    %angular noise difference
    angluar_diff = theta3-theta3_prime;
    %the negative loglikelihood of angle
    neg_ll_angle = log(2*pi) + log(besseli(0,nu_scaled)) - nu_scaled*cos(angluar_diff);

    %distance noise difference
    l3_prime = G3*h;
    dist_diff = l3-l3_prime;
    %the negative loglikelihood of distance on all trials 
    neg_ll_dist = 1/2*log(2*pi) + log(sigma_scaled) + (dist_diff^2)/(2*sigma_scaled^2);

    %total negative loglikelihood
    neg_ll = neg_ll_angle + neg_ll_dist;

    negloglikelihood = negloglikelihood + neg_ll;
end
end



