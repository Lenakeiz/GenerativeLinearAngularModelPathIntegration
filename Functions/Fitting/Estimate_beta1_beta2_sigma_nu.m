function [negloglikelihood] = Estimate_beta1_beta2_sigma_nu(beta1, beta2, sigma, nu, Input, config)
%   find the likelihood of the beta1 -beta2 - g2 - g3 - sigma - nu Model
%   Args:
%       beta1 is the decay factor for the mental distance on leg1
%       beta2 is the decay factor for the mental distance on leg2
%       g2 is the rotation gain for the second turn (measuring encoding error)
%       g3 is the rotation gain for the return (measuring production error)
%       sigma is the standard deviation for the Gaussian distribution of the return distance
%       nu is the standard deviation for the Gaussian distribution of the return angle
%       Input contains all the data information for estimating, see PerformGroupFit for how it was generated
%       config: self-explained

%% information necessary for running parameter estimation
DX              =   Input.DX;
THETAX          =   Input.THETADX;
flagOoB         =   Input.flagOoB;

sampleSize          =   size(DX,2);
negloglikelihood    =   0;

%% find the correct mean return angle based on all trials 
Alphas = zeros(sampleSize,1);
ActualAlphas = zeros(sampleSize,1);
Betas = zeros(sampleSize,1);
for tr = 1:sampleSize
    %extract the physical data info
    l1      = DX{tr}(1);
    l2      = DX{tr}(2);
    theta2  = THETAX{tr}(2);
    Betas(tr) = theta2;

    %calculate the correct return angle
    phy_p1  = [l1,0];
    phy_p2  = [l1+l2*cos(theta2),l2*sin(theta2)];
    vec1    = phy_p2-phy_p1; vec2 = [0,0]-phy_p2;
    alpha   = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
    alpha   = deg2rad(alpha);%transfer from degree to radians
    alpha   = mod(alpha, 2*pi);  %wrap to (0,2pi)  
    Alphas(tr) = alpha;
    
    %calculate the actuall return angle
    ActualAlphas(tr) = THETAX{tr}(3);
end
mean_angle = mean(ActualAlphas);

for tr = 1:sampleSize
    %% extract the physical data info
    l1          =       DX{tr}(1);
    l2          =       DX{tr}(2);
    l3          =       DX{tr}(3);
    theta2      =       THETAX{tr}(2); 
    theta3      =       THETAX{tr}(3); 
    
    %% whether to use weber's law to scaling the noise strength
    if config.useweber == true
        scale = sqrt(l1^2+l2^2); %noise scale
    else
        scale = 1; %noise scale
    end
    sigma_scaled = scale*sigma;
    nu_scaled = scale*nu;

    %mental point 1 (asuming a constant speed)
    %considering standing duration or not
    men_length1 = l1*beta1;
    men_p1 = [men_length1,0];
    
    theta2_prime = theta2;
    %theta2_prime = g2*theta2+mean_ecd_angle*(1-g2);

    %mental point 2, (asuming a constant speed)
    men_length2 = l2*beta2;
    men_p2      = [men_length1+men_length2*cos(theta2_prime),men_length2*sin(theta2_prime)];

    %calculate length of mental vector 3
    h           = norm(men_p2);
    %calculate turn angle of mental vector 3
    vec1        = men_p2-men_p1; 
    vec2        = [0,0]-men_p2;
    alpha       = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
    alpha       = deg2rad(alpha);   %transfer from degree to radians
    
    %mental turning angle
    theta3_prime = alpha;
    
    %angular noise difference
    angluar_diff = theta3-theta3_prime;

    neg_ll_angle = 1/2*log(2*pi) + log(nu_scaled) + (angluar_diff^2)/(2*nu_scaled^2); %Gaussian distribution

    %distance noise difference
    l3_prime    = h;
    dist_diff   = l3-l3_prime;

    %     %the negative loglikelihood of distance on non-OoB trials
    if flagOoB(tr)==0
        %this is a non-OoB trial
        neg_ll_dist = 1/2*log(2*pi) + log(sigma_scaled) + (dist_diff^2)/(2*sigma_scaled^2);
    else
        %this is an OoB trial
        neg_ll_dist = 0;
    end

    %total negative loglikelihood
    neg_ll = neg_ll_angle + neg_ll_dist;

    negloglikelihood = negloglikelihood + neg_ll;

end
end