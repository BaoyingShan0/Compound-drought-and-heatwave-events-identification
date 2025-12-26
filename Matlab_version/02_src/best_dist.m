function [ind_best, models]=best_dist(data)

% 2023.04.03
% best distribution for SPI, SHI, not for removing and merging

% aim: find the best distribution for data based on minimum AIC

% input:
% data: a sequence of variable

% output:
% ind_best: the index of the best-fitting in models
% models: all the distribution for chosen



%% set models for choosen
models=[ "Normal","Exponential","Gamma", "gev","InverseGaussian",...
    "logistic", "Loglogistic","Lognormal","Burr"];
k=[2, 1, 2, 3, 2, ...
    2, 2, 2, 3];% the number of parameter

%%
N=size(models, 2);

AIC=nan(N,1 ); 
for i=1:N

    try % if an error, then skip
        lastwarn('', ''); % reset the lastwarn message and id
        pd=fitdist(data, models(i) ) ;  % [shape, scale, location]
        % now if a warning was raised, warnMsg and warnId will not be empty.
        [~, warnId] = lastwarn();
        
        % you can check the warning message or id, or just throw the warning as an error if desired
        if(isempty(warnId)) % if there is no warning for fitting, then calculate MLE
            MLE=-negloglik(pd) ; %log likelihood
            AIC(i)=2*k(i)-2*MLE;            
        end
        
    catch
    end
end

% for return
[~,ind_best] = nanmin(AIC);



