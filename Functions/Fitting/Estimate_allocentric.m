function [negloglikelihood] = Estimate_allocentric(sigma, DX, THETAX, X, useweber)
%Estimate_ZL_mul means there exists a gain factor in the angle, which can not derive as a rotation matrix
%ESTIMATE Summary of this function goes here
%   sigma is the standard deviation for the Gaussian distribution of the return point
%   DX is a cell structure containing the segment of each trial
%   THETAX is the turning angle (wrong at the moment)
%   X is the data points

sampleSize = size(X,2);

negloglikelihood = 0;

for tr = 1:sampleSize
    if useweber
        l1 = DX{tr}(1); l2 = DX{tr}(2);
        scale = sqrt(l1^2+l2^2);
    else
        scale = 1;
    end
    %extract the physical end point from the data
    Phy_Endpoint = X{tr}(4,:); Phy_Endpoint=Phy_Endpoint';
    
    %calculate the negative log-likelihood center at 0. 
    % Note that the only parameter is sigma_prime.
    sigma_prime = scale*sigma;
    neg_ll = 1/2*log(2*pi) + log(sigma_prime) + sum(Phy_Endpoint.^2)/(2*sigma_prime^2);

    negloglikelihood = negloglikelihood + neg_ll;

end
end



