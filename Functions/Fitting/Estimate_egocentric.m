function [negloglikelihood] = Estimate_egocentric(sigma, nu, DX, THETAX, X, useweber)
%   sigma is the standard deviation for the Gaussian distribution of the return point
%   nu is the xxx
%   DX is a cell structure containing the segment of each trial
%   THETAX is the turning angle (wrong at the moment)
%   X is the data points
%   flagOoB is the flag of Out-of-Boundary trials, with 0 non-OoB trials 1 those OoB trials

sampleSize = size(X,2);

negloglikelihood = 0;

for tr = 1:sampleSize
    
    l1 = DX{tr}(1); l2 = DX{tr}(2); l3 = DX{tr}(3);

    theta2 = THETAX{tr}(2); theta3 = THETAX{tr}(3); 
    
    %physical endpoint 1
    phy_p1 = [l1,0];

    %physical endpoint 2
    phy_p2 = [l1+l2*cos(theta2), l2*sin(theta2)];

    %calculate length of physical vector 3
    h = norm(phy_p2);

    %calculate correct return angle 
    vec1 = phy_p2-phy_p1; vec2 = [0,0]-phy_p2;
    alpha = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
    %transfer from degree to radians
    alpha = deg2rad(alpha);
    
    %angular noise difference
    angluar_diff = theta3-alpha;
    %the negative loglikelihood of angle
    neg_ll_angle = log(2*pi) + log(besseli(0,nu)) - nu*cos(angluar_diff);

    %distance noise difference
    dist_diff = l3-h;
    %the negative loglikelihood of distance on all trials include OoB
    %trials
    neg_ll_dist = 1/2*log(2*pi) + log(sigma) + (dist_diff^2)/(2*sigma^2);
    
    %joint distribution (negative loglikelihood)
    neg_ll = neg_ll_angle + neg_ll_dist;

    negloglikelihood = negloglikelihood + neg_ll;

end
end



