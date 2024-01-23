clear; close all ;

p = 1;
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

%% plot lambda
% figure;
% for ir = 1:3
%     for con = [1 3 2 4]
%         subplot(2,2,con); hold on
%         histogram(fit.LambdaMinMSE{1,ir}(:,con), 0:10:100,'FaceAlpha',0.5);
%         xlabel('Lambda MinMSE');
%     end
% end
% figure;
% for ir = 1:3
%     for con = [1 3 2 4]
%         subplot(2,2,con); hold on
%         histogram(fit.Lambda1SE{1,ir}(:,con), 0:10:500,'FaceAlpha',0.5);
%         xlabel('Lambda MinMSE');
%     end
% end

%% plot individual beta
for ir = 1:6
    ic = 0;
for con = [1]
    figure;
    nfeatures = size(b_opt{1,con},2); ic = ic + 1;
    data = fit.Beta{1,ir}(1:nfeatures,:,con);
    regopt = fit_opt.Beta{1,ir}(1:nfeatures,:,con);
    
    for subi = 1:17
    subplot(3,6,subi); hold on
    
    if subi == 1  ylabel(strcat('Regularized beta by ', regressors(ir)));  end
    if con == 1 || con == 2
        xticks(xtk1); xlim([-0.5 nfeatures+1.5]); xticklabels(xt);
        ylim([-0.5 0.25]); yticks(-1.5:0.25:1.5);
    else
        xticks(xtk2); xlim([0.75 nfeatures+0.25]); xticklabels(xt);
        ylim([-0.85 0.4]); yticks(-0.8:0.2:0.4);
    end
    lineplot(0,'h','k--');
    
    p1 = plot(1:nfeatures, b_opt{1,con}, '--', 'Color','k','MarkerFaceColor','k','linewidth', 1);
    p1 = plot(1:nfeatures, regopt(:,subi), '-s', 'Color','k','MarkerFaceColor','k','linewidth', 1);    
    p4 = plot(1:nfeatures, data(:,subi), 's:', 'color', C22, 'MarkerFaceColor', C22,'linewidth', 1.5, 'MarkerSize', ms2);
    end
end
end

fa
%% plot beta
for ir = 1:7
    figure; ic = 0;
    for con = [1 3 2 4]
        nfeatures = size(b_opt{1,con},2); ic = ic + 1;
        
        subplot(2,2,Con(ic)); hold on
        data = fit.Beta{1,ir}(1:nfeatures,:,con);
        regopt = fit_opt.Beta{1,ir}(1:nfeatures,:,con);
        avbetaobs = mean(data,2); sebetaobs = std(data,[],2); % / sqrt(nsub) ;
        p1 = plot(1:nfeatures, b_opt{1,con}, '-', 'Color','k','linewidth', 1);
        p2 = plot(1:nfeatures, mean(regopt,2), '--s', 'Color','k','MarkerFaceColor','k','linewidth', 1);            
%         p3 = plot(1:nfeatures, data, '-', 'Color',C2);
        p3 = plot(1:nfeatures, data, 'wo', 'MarkerFaceColor',C2, 'MarkerSize',5);        
        p4 = errorbar(1:nfeatures, avbetaobs, sebetaobs, 's:', 'color', C22, 'MarkerFaceColor', C22,'linewidth', 1.5, 'MarkerSize', ms2, 'CapSize',cs);
        
        ylim([-1 0.75]); yticks(-1.5:0.5:1.5); title(tn(ic)); myfig;
        if con == 2  ylabel(strcat('Regularized beta by ', regressors(ir)));  end
        if con == 1 || con == 2
            xticks(xtk1); xlim([-0.5 nfeatures+1.5]); xticklabels(xt);
        else
            xticks(xtk2); xlim([0.75 nfeatures+0.25]); xticklabels(xt);
        end
    end
    legend([p1 p2 p4], {'Simulation','Special regression to S_i_d_e_a_l','Special regression to S_s_u_b'});
end

% average
mycolor = jet(7);
for i = 1:2                
figure; 
for ir = 1:6
    ic = 0;
    for con = [1 3 2 4]
        nfeatures = size(b_opt{1,con},2); ic = ic + 1;        
        if i == 1
            data = fit.Beta{1,ir}(1:nfeatures,:,con);
        else
            data = fit_opt.Beta{1,ir}(1:nfeatures,:,con);
        end
        
        subplot(2,2,Con(ic)); hold on
        avbetaobs = mean(data,2); sebetaobs = std(data,[],2); % / sqrt(nsub) ;
        p1 = plot(1:nfeatures, b_opt{1,con}, '-', 'Color','k');
        p2(ir,1) = plot(1:nfeatures, avbetaobs, 's:', 'color', mycolor(ir,:), 'MarkerEdgeColor', mycolor(ir,:),'MarkerFaceColor', mycolor(ir,:),'linewidth', 1.5, 'MarkerSize', ms2);
        
        ylim([-0.75 0.5]); yticks(-1.5:0.5:1.5); title(tn(ic));
        if con == 2  ylabel(strcat('Regularized beta by ', regressors(ir)));  end
        if con == 1 || con == 2
            xticks(xtk1); xlim([-0.5 nfeatures+1.5]); xticklabels(xt);
        else
            xticks(xtk2); xlim([0.75 nfeatures+0.25]); xticklabels(xt);
        end
    end
end
legend(p2,regressors);
end

%% ratio
for ir = 1:6
    figure; ic = 0;
    for con = [1 3 2 4]
        nfeatures = size(b_opt{1,con},2); ic = ic + 1;
        data = fit.Beta{1,ir}(1:nfeatures,:,con);
        
        betaratio = data ./ b_opt{1,con}' ;
        betaratio(betaratio==-inf) = NaN; betaratio(betaratio==inf) = NaN;
        isnanratio{1,ir}(1,con) = mean(mean(isnan(betaratio)));
        iszeroratio{1,ir}(1,con) = mean(mean(abs(betaratio)<=0.0001));         
        % Betaratio{1,ir}(:,:,con) = betaratio;
        avbetaratio = nanmean(betaratio,2); sebetaratio = nanstd(betaratio,[],2); % / sqrt(nsub) ;
        
        subplot(2,2,Con(ic)); hold on
        if con == 1 || con == 2
            xticks(xtk1); xlim([-0.5 nfeatures+1.5]); xticklabels(xt);
        else
            xticks(xtk2); xlim([0.75 nfeatures+0.25]); xticklabels(xt);
        end
        ylim([-20 20]); yticks(-20:5:20); title(tn(ic)); myfig;        
        lineplot(0, 'h','k--'); lineplot(1,'h', 'k--');

        p1 = plot(1:nfeatures, betaratio, 'wo', 'MarkerFaceColor',C2, 'MarkerSize',5);
        p2 = errorbar(1:nfeatures, avbetaratio, sebetaratio, 's:', 'color', C22, 'MarkerFaceColor', C22,'linewidth', 1.5, 'MarkerSize', ms2, 'CapSize',cs);
        
        ylabel('Beta_s_u_b / Beta_o_p_t');
    end
        legend(p2, regressors(ir));    
end



