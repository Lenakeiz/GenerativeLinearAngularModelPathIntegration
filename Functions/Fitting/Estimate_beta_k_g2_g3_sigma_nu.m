function [negloglikelihood] = Estimate_beta_k_g2_g3_sigma_nu(beta, k, g2, g3, sigma, nu, Input, config)
% Zilong Ji, UCL, zilong.ji@ucl.ac.uk
% Find the likelihood of the beta - k - g2 - g3 - sigma - nu Model
% Args:
% beta is the decay factor for the mental distance
% k is the velocity gain factor for the leaky integrator
% g2 is the rotation gain for the second turn (measuring encoding error)
% g3 is the rotation gain for the return (measuring production error)
% sigma is the standard deviation for the Gaussian distribution of the return distance
% nu is the standard deviation for the Gaussian distribution of the return angle
% Input contains all the data information for estimating, see PerformGroupFit for how it was generated
% For schematics about the steps of the generative linear and angular
% model please refer to Fig. 2 and online methods
% ===================================================================================

% Information necessary for running parameter estimation
DX              =   Input.DX;
THETAX          =   Input.THETADX;
L1Dur           =   Input.L1Dur;
L2Dur           =   Input.L2Dur;
StandingDur     =   Input.StandingDur;
flagOoB         =   Input.flagOoB;

sampleSize          =   size(DX,2);
negloglikelihood    =   0;

% Find the correct mean return angle and mean return distance
Alphas = zeros(sampleSize,1);
ActualAlphas = zeros(sampleSize,1);
Betas = zeros(sampleSize,1);
for tr = 1:sampleSize
    
    % Extract the real data info
    l1      = DX{tr}(1);
    l2      = DX{tr}(2);
    theta2  = THETAX{tr}(2);
    Betas(tr) = theta2;

    % Calculate the correct return angle
    phy_p1  = [l1,0];
    phy_p2  = [l1+l2*cos(theta2),l2*sin(theta2)];
    vec1    = phy_p2-phy_p1; vec2 = [0,0]-phy_p2;
    alpha   = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
    alpha   = deg2rad(alpha); %transfer from degree to radians
    alpha   = mod(alpha, 2*pi); %wrap to (0,2pi)  
    Alphas(tr) = alpha;
    
    % Calculate the actuall return angle
    ActualAlphas(tr) = THETAX{tr}(3);
end

mean_angle = mean(ActualAlphas);

for tr = 1:sampleSize
    % Extracting data from trial
    l1          =       DX{tr}(1);
    l2          =       DX{tr}(2);
    l3          =       DX{tr}(3);
    theta2      =       THETAX{tr}(2); 
    theta3      =       THETAX{tr}(3); 
    durationL1  =       L1Dur{tr}; 
    durationL2  =       L2Dur{tr};
    durationStand =     StandingDur{tr};
    
    % whether to use weber's law to scaling the noise strength
    if config.useweber == true
        scale = sqrt(l1^2+l2^2); %noise scale
    else
        scale = 1; %noise scale
    end
    sigma_scaled = scale*sigma;
    nu_scaled = scale*nu;

    % Calculation of mental point 1 (l1') 
    % Since we are using a leaky integration over time we considering
    % standing duration at cone 2 as well.
    if config.includeStand==true
        men_length1 = l1*k*(1-exp(-beta*durationL1))/(beta*durationL1)*exp(-beta*(durationL2+durationStand));
    else
        men_length1 = l1*k*(1-exp(-beta*durationL1))/(beta*durationL1)*exp(-beta*durationL2);
    end
    men_p1 = [men_length1,0];
    
    theta2_prime = g2*theta2;

    % Calculation of mental point 2
    men_length2 = l2*k*(1-exp(-beta*durationL2))/(beta*durationL2);
    men_p2      = [men_length1+men_length2*cos(theta2_prime),men_length2*sin(theta2_prime)];

    % Length of mental vector 3 (mental return to the origin)
    h           = norm(men_p2);
    
    % Turn angle of mental vector 3
    vec1        = men_p2-men_p1; 
    vec2        = [0,0]-men_p2;
    alpha       = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
    alpha       = deg2rad(alpha);   %transfer from degree to radians
    
    % Mental turning angle
    sign_alpha = sign(alpha);
    theta3_prime = g3*abs(alpha)+mean_angle*(1-g3); %reress to mean correct return angle
    theta3_prime = sign_alpha*theta3_prime;
    
    % Angular noise difference
    angluar_diff = theta3-theta3_prime;

    neg_ll_angle = 1/2*log(2*pi) + log(nu_scaled) + (angluar_diff^2)/(2*nu_scaled^2); %Gaussian distribution for angular return

    % Distance noise difference
    l3_prime    = h;
    dist_diff   = l3-l3_prime;

    % If OoB trial we set the contribution for the negative loglikelihood
    % of distance to zero
    if flagOoB(tr)==0
        % not OoB trial
        neg_ll_dist = 1/2*log(2*pi) + log(sigma_scaled) + (dist_diff^2)/(2*sigma_scaled^2);
    else
        % OoB trial
        neg_ll_dist = 0;
    end

    % Total negative loglikelihood
    neg_ll = neg_ll_angle + neg_ll_dist;

    negloglikelihood = negloglikelihood + neg_ll;

end

end