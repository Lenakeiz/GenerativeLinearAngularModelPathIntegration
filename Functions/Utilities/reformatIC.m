function [AIC, BIC, NLL] = reformatIC(IC)
%re-format IC from cell to array
%Args:
%   IC
%Returns:
%   AIC, BIC and Negative Loglikelihood
    numSubj = length(IC{1});
    AIC = zeros(numSubj,3); 
    BIC= zeros(numSubj,3); 
    NLL = zeros(numSubj,3);
    for cond = 1:3
        for subj = 1:numSubj
            AIC(subj, cond) = IC{cond}{subj}.aic;
            BIC(subj, cond) = IC{cond}{subj}.bic;
            NLL(subj, cond) = IC{cond}{subj}.negll;
        end
    end
    %remove nan rows from array
    AIC = removeNanRows(AIC);
    BIC = removeNanRows(BIC);
    NLL = removeNanRows(NLL);

    %average over conditions to get the IC mean
    AIC = mean(AIC,2);
    BIC = mean(BIC,2);
    NLL = mean(NLL,2);
end