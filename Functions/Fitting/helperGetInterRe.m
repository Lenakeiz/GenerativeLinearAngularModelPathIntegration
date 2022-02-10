function [l3_prime] = helperGetInterRe(G1, G2, G3, g2, g3, k3, sigma, nu, DX, THETAX, X)
%helper function to get intermediate results from the optimization target function
%   G1 is the gain of the length of the first leg 
%   G2 is the gain of the length of the second leg 
%   G3 is the gain of the length of the third leg 
%   g2 is the rotation angle from the direction of leg 2
%   g3 is the rotation angle from the direction of leg 3
%   k3 is the regression angle in leg 3
%   sigma is the standard deviation for the Gaussian distribution of the return point
%   nu is the xxx
%   DX is a cell structure containing the segment of each trial
%   THETAX is the turning angle (wrong at the moment)
%   X is the data points
%   flagOoB is the flag of Out-of-Boundary trials, with 0 non-OoB trials 1 those OoB trials

sampleSize = size(X,2);


for tr = 1:sampleSize
    
    l1 = DX{tr}(1); l2 = DX{tr}(2); l3 = DX{tr}(3);

    theta1 = THETAX{tr}(1); theta2 = THETAX{tr}(2); theta3 = THETAX{tr}(3); 

    %physical vector 3
    phy3 = [l3*cos(theta2+theta3),l3*sin(theta2+theta3)];

    %mental point 1
    men_p1 = [G1*l1,0];

    %mental point 2
    theta2_prime = g2*theta2;
    men_p2 = [G1*l1+G2*l2*cos(theta2_prime),G2*l2*sin(theta2_prime)];

    %calculate length of mental vector 3
    h = norm(men_p2);

    %distance noise difference
    l3_prime(tr) = G3*h;

    %calculate turn angle of mental vector 3
    vec1 = men_p2-men_p1; vec2 = [0,0]-men_p2;
    alpha = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
    %transfer from degree to radians
    alpha = deg2rad(alpha);
    %wrap to (0,2pi)
    alpha = mod(alpha, 2*pi);

    %considering execution turn error
    theta3_prime = g3*alpha+k3*(1-g3);
    %also wrap to (0,2pi)
    theta3_prime = mod(theta3_prime, 2*pi);
    
    %angular noise difference
    angluar_diff = theta3-theta3_prime;

end
end



