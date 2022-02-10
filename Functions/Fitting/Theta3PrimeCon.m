function [const, consteq] = Theta3PrimeCon(gamma, G3, g2, g3, b, DX, THETAX, X)
%%Nonlinear Constraints
%Nonlinear inequality constraints have the form c(x) ≤ 0, where c is a 
%vector of constraints, one component for each constraint. 
%Similarly, nonlinear equality constraints have the form ceq(x) = 0. 
%   gamma is the discount factor
%   G3 is the gain of the length of the third leg 
%   g2 is the rotation angle from the direction of leg 2
%   g3 is the rotation angle from the direction of leg 3
%   c is the systematic bias in the execution error
%   sigma is the standard deviation for the Gaussian distribution of the return point
%   nu is the xxx
%   DX is a cell structure containing the segment of each trial
%   THETAX is the turning angle (wrong at the moment)
%   X is the data points


sampleSize = size(X,2);

%get theta3_bar
theta3_bar_all = zeros(1,sampleSize);
for tr = 1:sampleSize
    theta3 = THETAX{tr}(3);
    theta3_bar_all(tr) = theta3;
end
theta3_bar = mean(theta3_bar_all);

for tr = 1:sampleSize
    
    l1 = DX{tr}(1); l2 = DX{tr}(2);

    theta2 = THETAX{tr}(2);

    %mental point 1
    %get G1
    G1 = gamma^2*G3;     
    men_p1 = [G1*l1,0];

    %mental point 2    
    theta2_prime = g2*theta2;
    %get G2
    G2 = gamma*G3;      
    men_p2 = [G1*l1+G2*l2*cos(theta2_prime),G2*l2*sin(theta2_prime)];

    %calculate turn angle of mental vector 3
    vec1 = men_p2-men_p1; vec2 = [0,0]-men_p2;
    alpha = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
    %transfer from degree to radians
    alpha = deg2rad(alpha);
    %wrap to (0,2pi)
    alpha = mod(alpha, 2*pi);

    %considering execution turn error
    theta3_prime = g3*alpha+b;

    const(2*tr) = -theta3_prime;
    const(2*tr+1) = theta3_prime-2*pi;

    consteq = [];

end
end