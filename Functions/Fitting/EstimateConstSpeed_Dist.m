function [NLL_Dist] = EstimateConstSpeed_Dist(beta, sigma, Input, config)
%   Distance parameter estimation
%   Args:
%       beta is the leaky integration decay factor
%       sigma is the standard deviation for the Gaussian distribution of the return point
%       Input contains all the data information for estimating, see PerformGroupFit for how it was generated
%       config: self-explained

%% information necessary for running parameter estimation
DX              =   Input.DX;
THETAX          =   Input.THETADX;
L1Dur           =   Input.L1Dur;
L2Dur           =   Input.L2Dur;
StandingDur     =   Input.StandingDur;
flagOoB         =   Input.flagOoB;

sampleSize      =   size(DX,2);
NLL_Dist        =   0;

for tr = 1:sampleSize
    %% extract the physical data info
    l1              =       DX{tr}(1);
    l2              =       DX{tr}(2);
    l3              =       DX{tr}(3);
    theta2          =       THETAX{tr}(2); 
    durationL1      =       L1Dur{tr}; 
    durationL2      =       L2Dur{tr};
    durationStand   =       StandingDur{tr};
    
    %% whether to use weber's law to scaling the noise strength
    if config.useweber == true
        scale = sqrt(l1^2+l2^2); %noise scale
    else
        scale = 1; %noise scale
    end
    sigma_scaled = scale*sigma;

    %mental point 1 (asuming a constant speed)

    %considering standing duration or not
    if config.includeStand==true
        men_length1 = l1*(1-exp(-beta*durationL1))/(beta*durationL1)*exp(-beta*(durationL2+durationStand));
    else
        men_length1 = l1*(1-exp(-beta*durationL1))/(beta*durationL1)*exp(-beta*durationL2);
    end
    men_p1 = [men_length1,0];
    
    theta2_prime = theta2; %considering no encoding error in the turning

    %mental point 2, (asuming a constant speed)
    men_length2 = l2*(1-exp(-beta*durationL2))/(beta*durationL2);
    men_p2      = [men_length1+men_length2*cos(theta2_prime),men_length2*sin(theta2_prime)];

    %calculate length of mental vector 3
    l3_prime    = norm(men_p2);

    %distance noise difference
    dist_diff   = l3-l3_prime;

    if config.useOoBTrial == true 
        %use OoB trials with the walking length replaced by the mean walking lenghth of good trials
        %the negative loglikelihood of distance on all trials 
        neg_ll_dist = 1/2*log(2*pi) + log(sigma_scaled) + (dist_diff^2)/(2*sigma_scaled^2);
    else
        %     %the negative loglikelihood of distance on non-OoB trials
        if flagOoB(tr)==0
            %this is a non-OoB trial
            neg_ll_dist = 1/2*log(2*pi) + log(sigma) + (dist_diff^2)/(2*sigma^2);
        else
            %this is an OoB trial
            neg_ll_dist = 0;
        end
    end
    
    %cumulating negative loglikelihood
    NLL_Dist = NLL_Dist + neg_ll_dist;
end
end



