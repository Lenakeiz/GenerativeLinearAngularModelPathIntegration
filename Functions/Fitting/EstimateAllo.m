function [negloglikelihood] = EstimateAllo(sigma, DX,  X, useweber)
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

    if useweber
        l1 = DX{tr}(1); l2 = DX{tr}(2);
        scale = sqrt(l1^2+l2^2); %noise scale
    else
        scale = 1; %noise scale
    end

    %the simplest generative model in an allcentric way
    %extract the physical end point from the data
    Phy_Endpoint = X{tr}(4,:); Phy_Endpoint=Phy_Endpoint';
    %calculate the negative log-likelihood center at 0. 
    % Note that the only parameter is sigma_prime.
    sigma_scaled = scale*sigma;
    neg_ll = 1/2*log(2*pi) + log(sigma_scaled) + sum(Phy_Endpoint.^2)/(2*sigma_scaled^2);                

    negloglikelihood = negloglikelihood + neg_ll;
end
end



