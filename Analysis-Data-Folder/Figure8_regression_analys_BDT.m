clear ; close all;

% Copy right: 12-3-2024. Keiji Ota (keiji.ota@nyu.edu)
% This is an incomplete code. It needs a large data file to run. 
% However, the idea is that we compute the optimal set point given a sample of N
% points. We then induce a subtle shift to a point in the sample and compute
% the sample variance induced by that shift. Then we recompute the optimal
% set point. We take a partial derivative of the optimal set point with
% respect to the induced shift in the sample. This numerical computation
% provides a regression coefficient of each sample point on the optimal set point. 
% Let me know if you need a complete code.

load egfunc; 

rng('shuffle');
nrep = 1000;
s = nan(nrep,4); 

tic
for con = 1:4
    if con == 1
        N = 30; eg = eg1;
    elseif con == 2
        N = 30; eg = eg2;
    elseif con == 3
        N = 5; eg = eg1;
    else
        N = 5; eg = eg2;
    end
    
    W = nan(nrep,N); x = [];   
    parfor irep = 1:nrep
        sigma = rand*10+10;
        X = randn(1,N)*sigma;
        X = sort(X,2) ;
        x(irep,:) = X;
        
        samplesigma = std(X,[],2);
        [S] = cal_Sstar(samplesigma,N,eg,prior_var,v_hut,Y);
        s(irep,con) = S;
        
        dXi = 1; shiftS = nan(1,N);
        for ix = 1:N
            XX = X; XX(1,ix) = X(1,ix) + dXi ;
            shiftsamplesigma = std(XX,[],2) ;
            [S] = cal_Sstar(shiftsamplesigma,N,eg,prior_var,v_hut,Y);
            shiftS(1,ix) = S ;
        end
        
        dS_dXi = shiftS - s(irep,con) ; % dS/dXi = dXi*Wi => Wi = dS/dXi / dXi
        W(irep,:) = dS_dXi / dXi ;
    end
    xx{1,con} = x;
    % subplot(2,2,con); plot(1:N, mean(W),'k-');
    meanW{1,con} = mean(W);
end
toc

save regression_results_BDT meanW xx s

