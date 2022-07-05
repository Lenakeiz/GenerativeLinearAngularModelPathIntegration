function [AIC, BIC, NLL] = reformatIC(IC, cond)
%Reformatting IC from cell to array for later figure plotting
%Args:
%   IC
%Returns:
%   AIC, BIC and Negative Loglikelihood
    
    if cond=="no change"
        index = 1;
    elseif cond=="no distal cue"
        index = 2;
    elseif cond=="no optical flow"
        index = 3;
    else
        %do nothing
    end

    numSubj = length(IC{1});
    AIC = zeros(numSubj,3); 
    BIC= zeros(numSubj,3); 
    NLL = zeros(numSubj,3);
    for cond_idx = 1:3
        for subj = 1:numSubj
            AIC(subj, cond_idx) = IC{cond_idx}{subj}.aic;
            BIC(subj, cond_idx) = IC{cond_idx}{subj}.bic;
            NLL(subj, cond_idx) = IC{cond_idx}{subj}.negll;
        end
    end
    %remove nan rows from array
    AIC = removeNanRows(AIC);
    BIC = removeNanRows(BIC);
    NLL = removeNanRows(NLL);
    
    if cond=="all"
        %average over conditions to get the IC mean
        AIC = mean(AIC,2);
        BIC = mean(BIC,2);
        NLL = mean(NLL,2);
    else
        AIC = AIC(:,index);
        BIC = BIC(:,index);
        NLL = NLL(:,index);   
    end
    
end