function [negloglikelihood] = Estimate(gamma, G3, g2, g3, b, sigma, nu, DX, THETAX, X, flagOoB)
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
%   flagOoB is the flag of Out-of-Boundary trials, with 0 non-OoB trials 1 those OoB trials

sampleSize = size(X,2);

negloglikelihood = 0;

%get theta3_bar
theta3_bar_all = zeros(1,sampleSize);
for tr = 1:sampleSize
    theta3 = THETAX{tr}(3);
    theta3_bar_all(tr) = theta3;
end
theta3_bar = mean(theta3_bar_all);

for tr = 1:sampleSize
    
    l1 = DX{tr}(1); l2 = DX{tr}(2); l3 = DX{tr}(3);

    theta1 = THETAX{tr}(1); theta2 = THETAX{tr}(2); theta3 = THETAX{tr}(3); 

    %physical vector 3
    phy3 = [l3*cos(theta2+theta3),l3*sin(theta2+theta3)];

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
    %also wrap to (0,2pi)
    %theta3_prime = mod(theta3_prime, 2*pi);
    
    %angular noise difference
    angluar_diff = theta3-theta3_prime;

    %the negative loglikelihood of angle
    neg_ll_angle = log(2*pi) + log(besseli(0,nu)) - nu*cos(angluar_diff);

    %distance noise difference
    l3_prime = G3*h;
    dist_diff = l3-l3_prime;

    %the negative loglikelihood of distance on all trials include OoB
    %trials
    neg_ll_dist = 1/2*log(2*pi) + log(sigma) + (dist_diff^2)/(2*sigma^2);

%     %the negative loglikelihood of distance on non-OoB trials
%     if flagOoB(tr)==0
%         %this is a non-OoB trial
%         neg_ll_dist = 1/2*log(2*pi) + log(sigma) + (dist_diff^2)/(2*sigma^2);
%     else
%         %this is an OoB trial
%         neg_ll_dist = 0;
%     end
    
    neg_ll = neg_ll_angle + neg_ll_dist;

    %neg_ll = 3/2*log(2*pi) + log(sigma) + log(besseli(0,nu)) + (dist_diff^2)/(2*sigma^2) -nu*cos(angluar_diff);

    %disp(num2str(tr)+" "+num2str(nu)+" "+num2str(sigma));

    negloglikelihood = negloglikelihood + neg_ll;

%     % merging angle and distance as below by a Gaussian noise
%     %execution mental vector 3. 
%     exe_men3 = [G3*h*cos(theta2+theta3_prime), G3*h*sin(theta2+theta3_prime)];
% 
%     %we want the execution mental vector 3 best explains the physica vector 3 
%     %by considering unseen factors as Gaussian noise 
%     %distribution, then we have negative loglikelihood of current trial
%     diff = phy3-exe_men3;
%     neg_ll = 1/2*log(2*pi) + log(sigma) + (sum(diff.^2))/(2*sigma^2);
%     negloglikelihood = negloglikelihood + neg_ll;
end
end



