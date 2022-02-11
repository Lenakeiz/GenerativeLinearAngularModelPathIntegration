function [negloglikelihood] = Estimate(gamma, G3, g2, g3, b, sigma, nu, DX, THETAX, X, ifallo, useweber)
%Estimate_ZL_mul means there exists a gain factor in the angle, which can not derive as a rotation matrix
%ESTIMATE Summary of this function goes here
%   gamma is the discount factor
%   G3 is the gain of the length of the third leg 
%   g2 is the rotation angle from the direction of leg 2
%   g3 is the rotation angle from the direction of leg 3
%   b is the systematic bias in the execution error
%   sigma is the standard deviation for the Gaussian distribution of the return point
%   nu is the xxx
%   DX is a cell structure containing the segment of each trial
%   THETAX is the turning angle (wrong at the moment)
%   X is the data points
%   ifallo 
%   useweber

sampleSize = size(X,2);

negloglikelihood = 0;

for tr = 1:sampleSize
    %extract the physical data info
    l1 = DX{tr}(1); l2 = DX{tr}(2); l3 = DX{tr}(3);
    theta2 = THETAX{tr}(2); theta3 = THETAX{tr}(3); 

    if useweber
        l1 = DX{tr}(1); l2 = DX{tr}(2);
        scale = sqrt(l1^2+l2^2); %noise scale
    else
        scale = 1; %noise scale
    end
    sigma_scaled = scale*sigma;
    nu_scaled = scale*nu;

    if ifallo ==false %normal way of calculating the likelihood
        %mental point 1
        %get G1
        G1 = gamma^2*G3; 
        men_p1 = [G1*l1,0];
    
        %mental point 2
        theta2_prime = g2*theta2;
        %get G2
        G2 = gamma*G3;    
        men_p2 = [G1*l1+G2*l2*cos(theta2_prime),G2*l2*sin(theta2_prime)];
    
        %calculate length of mental vector 3
        h = norm(men_p2);
        %calculate turn angle of mental vector 3
        vec1 = men_p2-men_p1; vec2 = [0,0]-men_p2;
        alpha = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
        %transfer from degree to radians
        alpha = deg2rad(alpha);
        %wrap to (0,2pi)
        alpha = mod(alpha, 2*pi);
    
        %considering execution turn error
        theta3_prime = g3*alpha+b;
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

    else %the simplest generative model in an allcentric way
        %extract the physical end point from the data
        Phy_Endpoint = X{tr}(4,:); Phy_Endpoint=Phy_Endpoint';
        %calculate the negative log-likelihood center at 0. 
        % Note that the only parameter is sigma_prime.
        sigma_scaled = scale*sigma;
        neg_ll = 1/2*log(2*pi) + log(sigma_scaled) + sum(Phy_Endpoint.^2)/(2*sigma_scaled^2);                
    end
    negloglikelihood = negloglikelihood + neg_ll;
end
end



