function [negloglikelihood] = EstimateDegradedLI_MeanReturnAng(beta, G3, g2, g3, b, sigma, nu, ProjSpeedL1, ProjSpeedL2, DX, THETAX, useweber)
%   EstimateLI means using the degraded leaky integration model when estimating the parameters
%   degraded means a constant speed
%   ESTIMATELI Summary of this function goes here:
%   beta is the leaky integration decay factor
%   G3 is the gain of the length of the third leg 
%   g2 is the rotation angle from the direction of leg 2
%   g3 is the rotation angle from the direction of leg 3
%   sigma is the standard deviation for the Gaussian distribution of the return point
%   nu decribes the noise strength in the Von Mises distribution
%   ProjSpeedL1, speed projected onto the first outbound path
%   ProjSpeedL2, speed projected onto the second outbound path
%   DX is a cell structure containing the segment of each trial
%   THETAX is the turning angle
%   useweber decides whether the noise scales with the walking distance 

sampleSize = size(DX,2);
deltat = 0.1; %the recording interval, always 0.1s.
negloglikelihood = 0;

%find the correct mean return angle based on all trials 
Alphas = zeros(sampleSize,1);
for tr = 1:sampleSize
    %extract the physical data info
    l1 = DX{tr}(1);l2 = DX{tr}(2);
    theta2 = THETAX{tr}(2);

    %calculate the correct return angle
    phy_p1 = [l1,0];
    phy_p2 = [l1+l2*cos(theta2),l2*sin(theta2)];
    vec1 = phy_p2-phy_p1; vec2 = [0,0]-phy_p2;
    alpha = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
    %transfer from degree to radians
    alpha = deg2rad(alpha);
    %wrap to (0,2pi)
    alpha = mod(alpha, 2*pi);    
    Alphas(tr) = alpha;
end
mean_angle = mean(Alphas);

for tr = 1:sampleSize
    %extract the physical data info
    l1 = DX{tr}(1);l2 = DX{tr}(2);l3 = DX{tr}(3);
    theta2 = THETAX{tr}(2); theta3 = THETAX{tr}(3); 

    timeL1 = ProjSpeedL1{1, tr}; 
    timeIntervalL1 = [diff(timeL1);deltat]; %add 0.1 to keep the dimension the same
    durationL1 = sum(timeIntervalL1); %total duration of walking in leg 1
    timeCumL1 = cumsum(timeIntervalL1); %cumulative sum of the time interval, for later use
    speedL1 = ProjSpeedL1{2, tr};
    sumV1 = sum(speedL1*deltat); %replace l1 with sumV1, l1 is the perfect distance of outbound path1

    timeL2 = ProjSpeedL2{1, tr}; 
    timeIntervalL2 = [diff(timeL2);0.1]; %add 0.1 to keep the dimension the same
    durationL2 = sum(timeIntervalL2); %total duration of walking in leg 2
    timeCumL2 = cumsum(timeIntervalL2); %cumulative sum of the time interval, for later use
    speedL2 = ProjSpeedL2{2, tr};
    sumV2 = sum(speedL2*deltat); %replace l1 with sumV1, l1 is the perfect distance of outbound path1

    if useweber
        l1 = DX{tr}(1); l2 = DX{tr}(2);
        scale = sqrt(l1^2+l2^2); %noise scale
    else
        scale = 1; %noise scale
    end
    sigma_scaled = scale*sigma;
    nu_scaled = scale*nu;

    %mental point 1, asuming a constant speed, which will integrate upto l1 in duration1
    men_length1 = l1*(1-exp(-beta*durationL1))/(beta*durationL1)*exp(-beta*durationL2);

    men_p1 = [men_length1,0];
    
    theta2_prime = g2*theta2;

    %mental point 2, asuming a constant speed, which will integrate upto l2
    %in duration2
    men_length2 = l2*(1-exp(-beta*durationL2))/(beta*durationL2);
    men_p2 = [men_length1+men_length2*cos(theta2_prime),men_length2*sin(theta2_prime)];

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
    theta3_prime = g3*alpha+(1-g3)*mean_angle;
    %also wrap theta3_prime to (0,2pi) (We don't have to coz of the 'cos' below, 
    %also coz of the nonconstriant (if there is) which restrict theta3_prime in [0,2pi]).
    theta3_prime = mod(theta3_prime, 2*pi);
    
    %angular noise difference
    angluar_diff = theta3-theta3_prime;

    %the negative loglikelihood of angle
    neg_ll_angle = log(2*pi) + log(besseli(0,nu_scaled)) - nu_scaled*cos(angluar_diff);
    %neg_ll_angle = 1/2*log(2*pi) + log(nu_scaled) + (angluar_diff^2)/(2*nu_scaled^2);

    %distance noise difference
    l3_prime = G3*h;
    dist_diff = l3-l3_prime;

    %the negative loglikelihood of distance on all trials 
    neg_ll_dist = 1/2*log(2*pi) + log(sigma_scaled) + (dist_diff^2)/(2*sigma_scaled^2);

    %total negative loglikelihood
    neg_ll = neg_ll_angle + neg_ll_dist;

    negloglikelihood = negloglikelihood + neg_ll;
end
end



