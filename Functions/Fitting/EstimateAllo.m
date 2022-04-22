function [negloglikelihood] = EstimateAllo(sigma, Input, config)
%   The simplest generative model in an allcentric way
%   Args:
%       sigma is the standard deviation for the Gaussian distribution of the return point
%       Input contains all the data information for estimating, see PerformGroupFit for how it was generated
%       config: self-explained

DX  =  Input.DX;
X   =  Input.X;

sampleSize = size(X,2);
negloglikelihood = 0;

for tr = 1:sampleSize
    
    %% whether to use weber's law to scaling the noise strength
    if config.useweber==true
        l1 = DX{tr}(1); l2 = DX{tr}(2);
        scale = sqrt(l1^2+l2^2); %noise scale
    else
        scale = 1; %noise scale
    end

    %% extract the physical end point from the data
    Phy_Endpoint = X{tr}(4,:); Phy_Endpoint=Phy_Endpoint';
    %% calculate the negative log-likelihood center at 0. 
    % Note that the only parameter is sigma_prime.
    sigma_scaled = scale*sigma;
    neg_ll = 1/2*log(2*pi) + log(sigma_scaled) + sum(Phy_Endpoint.^2)/(2*sigma_scaled^2);                

    negloglikelihood = negloglikelihood + neg_ll;
end
end



