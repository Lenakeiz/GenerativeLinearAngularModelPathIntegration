function [NLL_Angle] = EstimateConstSpeed_AngleRGb(g3, b, nu, beta, Input, config)
%   Angle parameter estimation
%   Args:    
%       g3 is the rotation angle from the direction of leg 3
%       b is the systematic bias in the execution error
%       nu decribes the noise strength in the Von Mises distribution
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

sampleSize      =   size(DX,2);
NLL_Angle       =   0;

for tr = 1:sampleSize
    %% extract the physical data info
    l1          =       DX{tr}(1);
    l2          =       DX{tr}(2);
    theta2      =       THETAX{tr}(2); 
    theta3      =       THETAX{tr}(3); 
    durationL1  =       L1Dur{tr}; 
    durationL2  =       L2Dur{tr};
    durationStand =     StandingDur{tr};
    
    %% whether to use weber's law to scaling the noise strength
    if config.useweber == true
        scale = sqrt(l1^2+l2^2); %noise scale
    else
        scale = 1; %noise scale
    end
    nu_scaled = scale*nu;

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

    %calculate turn angle of mental vector 3
    vec1        = men_p2-men_p1; 
    vec2        = [0,0]-men_p2;
    alpha       = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
    alpha       = deg2rad(alpha);   %transfer from degree to radians
    alpha       = mod(alpha, 2*pi); %wrap to (0,2pi)

    %regress to the angle b which is also a parameter
    theta3_prime = g3*alpha+(1-g3)*b;

    %also wrap theta3_prime to (0,2pi) (We don't have to coz of the 'cos' below, 
    %also coz of the nonconstriant (if there is) which restrict theta3_prime in [0,2pi]).
    theta3_prime = mod(theta3_prime, 2*pi);
    
    %angular noise difference
    angluar_diff = theta3-theta3_prime;
    %the negative loglikelihood of angle
    neg_ll_angle = log(2*pi) + log(besseli(0,nu_scaled)) - nu_scaled*cos(angluar_diff);
    
    %cumulating negative loglikelihood
    NLL_Angle = NLL_Angle + neg_ll_angle;
end
end



