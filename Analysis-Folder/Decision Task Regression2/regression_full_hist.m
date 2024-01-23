clear; close all ;

p = 2;
if p == 1
    load regression_result_X1X30;  % observation with 500 bootsamples
elseif p == 2
    load regression_result2_X1X30; % observation with 100? bootsamples
end
xt = {'Pmin','P25%','P50%','P75%','PMax'}; xtk1 = [1, 8, 15, 23, 30]; xtk2 = 1:5;
load sim_weight;

b_opt = meanW; nsub = 17;
tn = {'30 dots, 0 penalty', '5 dots, 0 penalty', '30 dots, -500 penalty','5 dots, -500 penalty'};
C1 = [150 150 150]./255 ; C11 = [10 10 10]./255 ;
C2 = [255 200 200]./255 ; C22 = [255 75 75]./255 ;
ms = 8; ms2 = 5; cs = 10; Con = 1:4;

regressors = {'LASSO min MSE','LASSO 1SE','Elastic Net (Alpha = 0.01), min MSE',...,
    'Elastic Net (Alpha = 0.01), 1SE','Ridge min MSE','Ridge 1SE','LM'};
fit.Beta{1,1} = fit.BetaMinMSE{1,1}; fit.Beta{1,2} = fit.Beta1SE{1,1};
fit.Beta{1,3} = fit.BetaMinMSE{1,2}; fit.Beta{1,4} = fit.Beta1SE{1,2};
fit.Beta{1,5} = fit.BetaMinMSE{1,3}; fit.Beta{1,6} = fit.Beta1SE{1,3};
fit.Beta{1,7} = fit.BetaMinMSE{1,4};

fit_opt.Beta{1,1} = fit_opt.BetaMinMSE{1,1}; fit_opt.Beta{1,2} = fit_opt.Beta1SE{1,1};
fit_opt.Beta{1,3} = fit_opt.BetaMinMSE{1,2}; fit_opt.Beta{1,4} = fit_opt.Beta1SE{1,2};
fit_opt.Beta{1,5} = fit_opt.BetaMinMSE{1,3}; fit_opt.Beta{1,6} = fit_opt.Beta1SE{1,3};
fit_opt.Beta{1,7} = fit_opt.BetaMinMSE{1,4};


%% ratio
ir = 5;
for i = 1:6
figure; ic = 0;
for con = [1 3 2 4]
    nfeatures = size(b_opt{1,con},2); ic = ic + 1;  x = 1:nfeatures;
    data = fit.Beta{1,ir}(1:nfeatures,:,con);
    regopt = fit_opt.Beta{1,ir}(1:nfeatures,:,con);

    betaratio1 = data ./ b_opt{1,con}' ;
    betaratio2 = data ./ regopt; 
    Betaratio{1,con} = betaratio1;

    tmp = abs(b_opt{1,con}') >= 0.01;    
    
    subplot(2,2,Con(ic)); hold on
    if con == 1 || con == 2
        xticks(xtk1); xlim([-0.5 nfeatures+1.5]); xticklabels(xt);
    else
        xticks(xtk2); xlim([0.75 nfeatures+0.25]); xticklabels(xt);
    end
    ylim([-10 10]); yticks(-10:1:10); title(tn(ic)); myfig;
    lineplot(0, 'h','k--'); lineplot(1,'h', 'k--');
    ylabel('Beta_s_u_b / Beta_o_p_t');
        
    if i == 1
        avbetaratio = mean(betaratio1,2); sebetaratio = nanstd(betaratio1,[],2); % / sqrt(nsub) ;
        p2 = errorbar(x(tmp), avbetaratio(tmp,1), sebetaratio(tmp,1), 's:', 'color', C22, 'MarkerFaceColor', C22,'linewidth', 1.5, 'MarkerSize', ms2, 'CapSize',cs);
    elseif i == 2
        avbetaratio = mean(betaratio1,2); sebetaratio = 2*nanstd(betaratio1,[],2) / sqrt(nsub) ;
        p2 = errorbar(x(tmp), avbetaratio(tmp,1), sebetaratio(tmp,1), 's:', 'color', C22, 'MarkerFaceColor', C22,'linewidth', 1.5, 'MarkerSize', ms2, 'CapSize',cs);
    elseif i == 3        
        medianratio = median(betaratio1,2); sratio = sort(betaratio1,2); pctl25 = sratio(:,4); pctl75 = sratio(:,13);
        p2 = errorbar(x(tmp), medianratio(tmp,1), pctl25(tmp,1),pctl75(tmp,1), 's:', 'color', C22, 'MarkerFaceColor', C22,'linewidth', 1.5, 'MarkerSize', ms2, 'CapSize',cs);
    elseif i == 4
        avbetaratio = mean(betaratio2,2); sebetaratio = nanstd(betaratio2,[],2); % / sqrt(nsub) ;
        p2 = errorbar(x(tmp), avbetaratio(tmp,1), sebetaratio(tmp,1), 's:', 'color', C22, 'MarkerFaceColor', C22,'linewidth', 1.5, 'MarkerSize', ms2, 'CapSize',cs);
    elseif i == 5
        avbetaratio = mean(betaratio2,2); sebetaratio = 2*nanstd(betaratio2,[],2) / sqrt(nsub) ;
        p2 = errorbar(x(tmp), avbetaratio(tmp,1), sebetaratio(tmp,1), 's:', 'color', C22, 'MarkerFaceColor', C22,'linewidth', 1.5, 'MarkerSize', ms2, 'CapSize',cs);
    elseif i == 6        
        medianratio = median(betaratio2,2); sratio = sort(betaratio2,2); pctl25 = sratio(:,4); pctl75 = sratio(:,13);
        p2 = errorbar(x(tmp), medianratio(tmp,1), pctl25(tmp,1),pctl75(tmp,1), 's:', 'color', C22, 'MarkerFaceColor', C22,'linewidth', 1.5, 'MarkerSize', ms2, 'CapSize',cs);                
    end
end
end 

%% histogram
for con = [1 2 3 4]
    figure;
    if con <= 2
        for ip = 1:30
            data = Betaratio{1,con}(ip,:);
            subplot(5,6,ip)
            hist(data);
        end
    else
        for ip = 1:5
            data = Betaratio{1,con}(ip,:);
            subplot(1,5,ip)
            hist(data);
        end
    end
end

%% 95% CI
figure;
ic = 0; bootN=10000; subn = 17;
for con = 1
    nfeatures = size(b_opt{1,con},2); ic = ic + 1;  x = 1:nfeatures;
    data = fit.Beta{1,ir}(1:nfeatures,:,con);
    regopt = fit_opt.Beta{1,ir}(1:nfeatures,:,con);

    betaratio = data ./ b_opt{1,con}' ;
    
    meanratio = mean(betaratio,2); twoSE = 1.96*std(betaratio,[],2)/sqrt(subn);
    
    bootmean = [];
    for ib = 1:bootN    
        tmp = randi(subn,subn,1);
        bootmean(ib,:) = mean(betaratio(:,tmp),2) ;
    end
    sbootmean = sort(bootmean);    
    pctl025 = sbootmean(bootN*0.025,:); pctl975 = sbootmean(bootN*0.975,:);    
    for ix = 1:30
        subplot(5,6,ix)
        histogram(bootmean(:,ix),'Normalization','Probability'); hold on
        lineplot(meanratio(ix), 'v','k:'); plot([meanratio(ix)-twoSE(ix), meanratio(ix)+twoSE(ix)], ones(2,1)*0.05,'k-');
        lineplot(mean(bootmean(:,ix)), 'v', 'k--'); plot([pctl025(ix) pctl975(ix)], ones(2,1)*0.055,'r-');
        xlim([min(bootmean(:,ix)) max(bootmean(:,ix))]);
        if ix == 1    ylabel('Distribution of mean of ratio'); end
    end
end












