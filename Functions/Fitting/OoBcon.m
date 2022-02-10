function [c, ceq] = OoBcon(G1, G2, G3, g2, g3, k3, sigma, nu, DX, THETAX, X, OoBLen, flagOoB)
%%Nonlinear Constraints
%Nonlinear inequality constraints have the form c(x) ≤ 0, where c is a 
%vector of constraints, one component for each constraint. 
%Similarly, nonlinear equality constraints have the form ceq(x) = 0. 
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
%   OoBLen is the Out-of-Boundary length with 0 for non-oob trials, >0 for oob trials.
%   flagOoB is the flag of Out-of-Boundary trials, with 0 non-OoB trials 1 those OoB trials


sampleSize = size(X,2);


for tr = 1:sampleSize
    
    l1 = DX{tr}(1); l2 = DX{tr}(2); l3 = DX{tr}(3);

    theta2 = THETAX{tr}(2);

    %mental point 1
    men_p1 = [G1*l1,0];

    %mental point 2
    theta2_prime = g2*theta2;
    men_p2 = [G1*l1+G2*l2*cos(theta2_prime),G2*l2*sin(theta2_prime)];

    %calculate length of mental vector 3
    h = norm(men_p2);

    l3_prime = G3*h;

    %set the nonlinear constraint, can simply write as l3_hat=OoBLen(tr),
    %here write as if just for clear understanding
    if flagOoB(tr)==0
        %this is a non-OoB trial
        l3_hat = 0; %this constraintion will always meet
    else
        %this is an OoB trial
        l3_hat = OoBLen(tr); 
    end
    c(tr) = l3_hat-l3_prime;

    ceq = [];

end
end