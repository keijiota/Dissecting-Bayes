clear; close all ;

p = 2;
if p == 1
    load regression_result_X1X30;  % observation with 500 bootsamples
elseif p == 2
    load regression_result2_X1X30; % observation with 100? bootsamples
end
xt1 = {'Pmin','','','','','','','P25%','','','','','','','P50%',...,
       '','','','','','','','P75%','','','','','','','PMax'}; xtk1 = 1:30; % [1, 8, 15, 23, 30] 
xt11 = {'Pmin/Pmax','','','','','','','P25%/P75%','','','','','','','P50%'}; xtk11 = 1:15; 
xt2 = {'Pmin','P25%','P50%','P75%','PMax'}; xtk2 = 1:5;
xt22 = {'Pmin/PMax','P25%/75%','P50%'}; xtk22 = 1:3;

load sim_weight;

b_opt = meanW; nsub = 17;
tn  = {'30 dots, 0 penalty', '5 dots, 0 penalty', '30 dots, -500 penalty','5 dots, -500 penalty'};

C = [80 180 255;
     25 25 255;     
     255 140 0;      
     200 0 0]./255 ; 
C22 = [255 75 75]./255 ;

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


%% plot beta
ir = 5;
figure; ic = 0;
for con = [1 3 2 4]
    nfeatures = size(b_opt{1,con},2); ic = ic + 1;
    
    subplot(2,2,Con(ic)); hold on
    data = fit.Beta{1,ir}(1:nfeatures,:,con);
    regopt = fit_opt.Beta{1,ir}(1:nfeatures,:,con);
    avbetaobs = mean(data,2); sebetaobs = std(data,[],2) / sqrt(nsub) ;
    
    title(tn(ic)); myfig;
    if con == 2  ylabel(strcat('Regularized beta by ', regressors(ir)));  end
    if con == 1 || con == 2
        xticks(xtk1); xlim([-0.5 nfeatures+1.5]); xticklabels(xt1);
        ylim([-0.35 0.15]); yticks(-1.5:0.1:1.5);
    else
        xticks(xtk2); xlim([0.75 nfeatures+0.25]); xticklabels(xt2);
        ylim([-0.65 0.3]); yticks(-0.8:0.2:0.4);
    end
    lineplot(0,'h','k--');
    
    p1 = plot(1:nfeatures, b_opt{1,con}, '-s', 'Color','k','MarkerFaceColor','k','linewidth', 1);
%     p2 = plot(1:nfeatures, regopt, '-s', 'Color','k','MarkerFaceColor','k','linewidth', 1);
    %     p2 = plot(1:nfeatures, mean(regopt,2), '--s', 'Color','k','MarkerFaceColor','k','linewidth', 1);
    p4 = errorbar(1:nfeatures, avbetaobs, sebetaobs, 's:', 'color', C22, 'MarkerFaceColor', C22,'linewidth', 1.5, 'MarkerSize', ms2, 'CapSize',cs);
    
end
legend([p1 p4], {'Simulation of ideal weights','Subjects weights'});

%% ratio without error bars
figure; ic = 0; ymin = -3.5; ymax = 5.5;
for con = [1 3 2 4]
    nfeatures = size(b_opt{1,con},2); ic = ic + 1;  x = 1:nfeatures;
    data = fit.Beta{1,ir}(1:nfeatures,:,con);
    regopt = fit_opt.Beta{1,ir}(1:nfeatures,:,con);
    
    betaratio = data ./ b_opt{1,con}' ;
    %     betaratio2 = data ./ regopt;
    betaratio(betaratio==-inf) = NaN; betaratio(betaratio==inf) = NaN;
    isnanratio{1,ir}(1,con)  = mean(mean(isnan(betaratio)));
    iszeroratio{1,ir}(1,con) = mean(mean(abs(betaratio)<=0.0001));
    Betaratio{1,con} = betaratio;
    
    tmp = abs(b_opt{1,con}') >= 0.01; 
    
    if con == 1 || con == 2
        subplot(1,2,1); hold on        
        xticks(xtk1); xlim([-0.5 nfeatures+1.5]);  xticklabels(xt1);
        tmp = ones(1,30); tmp(13:17) = 0; % at least 13-17
    else
        subplot(1,2,2); hold on
        xticks(xtk2); xlim([0.75 nfeatures+0.25]); xticklabels(xt2);
        tmp = ones(1,5); tmp(3) = 0;         
    end
    ismask = find(tmp == 0);
    ylim([ymin ymax]); yticks(-10:1:10); myfig;
    lineplot(0, 'h','k--'); lineplot(1,'h', 'k--');
    ylabel('Beta_s_u_b / Beta_o_p_t');
    
    avbetaratio = mean(betaratio,2);     
    p1(con) = plot(1:min(ismask)-1, avbetaratio(1:min(ismask)-1,1), 's:', 'color', C(con,:), 'MarkerFaceColor', C(con,:),'linewidth', 1.5, 'MarkerSize', ms2);
    p2 = plot(max(ismask)+1:nfeatures, avbetaratio(max(ismask)+1:nfeatures,1), 's:', 'color', C(con,:), 'MarkerFaceColor', C(con,:),'linewidth', 1.5, 'MarkerSize', ms2);
    rectangle('Position',[min(ismask)-0.5,ymin,max(ismask)-min(ismask)+1,ymax-ymin],'FaceColor',[0.9 0.9 0.9],'EdgeColor','w');       
end
subplot(1,2,1); legend(p1([1 2]), tn([1 3]));
subplot(1,2,2); legend(p1([3 4]), tn([2 4]));

%% ratio with error bars
% figure; ic = 0;  ymin = -5; ymax = 8; P = [0 0.5 0 0.5];
% for con = [1 3 2 4]
%     nfeatures = size(b_opt{1,con},2); ic = ic + 1;  x = 1:nfeatures;    
%     betaratio = Betaratio{1,con};    
%     tmp = abs(b_opt{1,con}') >= 0.01;
%     
%     if con == 1 || con == 2
%         subplot(1,2,1); hold on        
%         xticks(xtk1); xlim([-0.5 nfeatures+1.5]);  xticklabels(xt1); p = P(con);
%         tmp = ones(1,30); tmp(13:17) = 0; % at least 13-17
%     else
%         subplot(1,2,2); hold on
%         xticks(xtk2); xlim([0.75 nfeatures+0.25]); xticklabels(xt2); p = P(con);
%         tmp = ones(1,5); tmp(3) = 0;  
%     end
%     ismask = find(tmp == 0);
%     ylim([ymin ymax]); yticks(-10:1:10); title(tn(ic)); myfig;
%     lineplot(0, 'h','k--'); lineplot(1,'h', 'k--');
%     ylabel('Beta_s_u_b / Beta_o_p_t');
%     
%     avbetaratio = mean(betaratio,2); sebetaratio = nanstd(betaratio,[],2) / sqrt(nsub) ;
%     p1(con) = errorbar(1+p:min(ismask)-1+p, avbetaratio(1:min(ismask)-1,1), sebetaratio(1:min(ismask)-1,1),'s:', 'color', C(con,:), 'MarkerFaceColor', C(con,:),'linewidth', 1.5, 'MarkerSize', ms2, 'CapSize',cs);
%     p2      = errorbar(max(ismask)+1+p:nfeatures+p, avbetaratio(max(ismask)+1:nfeatures,1), sebetaratio(max(ismask)+1:nfeatures,1), 's:', 'color', C(con,:), 'MarkerFaceColor', C(con,:),'linewidth', 1.5,'MarkerSize', ms2, 'CapSize',cs);
%     rectangle('Position',[min(ismask)-0.5,ymin,max(ismask)-min(ismask)+1,ymax-ymin],'FaceColor',[0.9 0.9 0.9],'EdgeColor','w');       
% end
% subplot(1,2,1); legend(p1([1 2]), tn([1 3]));
% subplot(1,2,2); legend(p1([3 4]), tn([2 4]));

figure; ic = 0;  ymin = -5; ymax = 7.5;
for con = [1 3 2 4]
    nfeatures = size(b_opt{1,con},2); ic = ic + 1;  x = 1:nfeatures;    
    betaratio = Betaratio{1,con};    
    tmp = abs(b_opt{1,con}') >= 0.01;
    
    subplot(2,2,Con(ic)); hold on
    if con == 1 || con == 2
        xticks(xtk1); xlim([-0.5 nfeatures+1.5]);  xticklabels(xt1);
        tmp = ones(1,30); tmp(14:17) = 0; % at least 13-17
    else
        xticks(xtk2); xlim([0.75 nfeatures+0.25]); xticklabels(xt2);
        tmp = ones(1,5); tmp(3) = 0;  
    end
    ismask = find(tmp == 0);
    ylim([ymin ymax]); yticks(-10:2.5:10); title(tn(ic)); myfig;
    lineplot(0, 'h','k--'); lineplot(1,'h', 'k--');
    ylabel('Beta_s_u_b / Beta_o_p_t');
    
    avbetaratio = mean(betaratio,2); sebetaratio = nanstd(betaratio,[],2) / sqrt(nsub) ;
    p1(con) = errorbar(1:min(ismask)-1, avbetaratio(1:min(ismask)-1,1), sebetaratio(1:min(ismask)-1,1),'s:', 'color', C(con,:), 'MarkerFaceColor', C(con,:),'linewidth', 1.5, 'MarkerSize', ms2, 'CapSize',cs);
    p2      = errorbar(max(ismask)+1:nfeatures, avbetaratio(max(ismask)+1:nfeatures,1), sebetaratio(max(ismask)+1:nfeatures,1), 's:', 'color', C(con,:), 'MarkerFaceColor', C(con,:),'linewidth', 1.5,'MarkerSize', ms2, 'CapSize',cs);    
    rectangle('Position',[min(ismask)-0.5,ymin,max(ismask)-min(ismask)+1,ymax-ymin],'FaceColor',[0.9 0.9 0.9],'EdgeColor','w');       
end

figure; ic = 0; 
for con = [1 3 2 4]
    nfeatures = size(b_opt{1,con},2); ic = ic + 1;  x = 1:nfeatures;    
    betaratio = Betaratio{1,con};    
    tmp = abs(b_opt{1,con}') >= 0.01;
    
    subplot(2,2,Con(ic)); hold on
    if con == 1 || con == 2
        xticks(xtk11+0.125); xlim([-0.5 15+1.5]);  xticklabels(xt11);
        tmp = ones(1,30); tmp(14:17) = 0; % at least 13-17
    else
        xticks(xtk22+0.125); xlim([0.75 3+0.25]); xticklabels(xt22);
        tmp = ones(1,5); tmp(3) = 0;  
    end
    ismask = find(tmp == 0);
    ylim([ymin ymax]); yticks(-10:2.5:10); title(tn(ic)); myfig;
    lineplot(0, 'h','k--'); lineplot(1,'h', 'k--');
    ylabel('Beta_s_u_b / Beta_o_p_t');
    
    avbetaratio = mean(betaratio,2); sebetaratio = nanstd(betaratio,[],2) / sqrt(nsub) ;
    p1(con) = errorbar(1:min(ismask)-1, avbetaratio(1:min(ismask)-1,1), sebetaratio(1:min(ismask)-1,1),'s:', 'color', C(con,:), 'MarkerFaceColor', C(con,:),'linewidth', 1.5, 'MarkerSize', ms2, 'CapSize',cs);
    p2      = errorbar(1+0.25:min(ismask)-1+0.25, flip(avbetaratio(max(ismask)+1:nfeatures,1)), flip(sebetaratio(max(ismask)+1:nfeatures,1)), 's:', 'color', C(con,:), 'MarkerFaceColor', C(con,:),'linewidth', 1.5,'MarkerSize', ms2, 'CapSize',cs);
    rectangle('Position',[min(ismask)-0.5,ymin,max(ismask)-min(ismask)+1,ymax-ymin],'FaceColor',[0.9 0.9 0.9],'EdgeColor','w');       
end

%% histogram
for con = [1 2 3 4]
    figure; i = 0;
    for ip = 1:size(Betaratio{1,con},1)
        data = Betaratio{1,con}(ip,:); d = mean(abs(diff(data)))-min(abs(diff(data)));
        subplot(5,6,ip)
        histogram(data, -2000:d:2000,'Normalization','Probability');
        if con == 1 || con == 2
            if ip == 1 || ip == 8 || ip == 15 || ip == 23 || ip == 30
                i = i+1; title(xt2(i)); xlabel('Beta_s_u_b / Beta_o_p_t');
            end
        else
            i = i+1; title(xt2(i)); xlabel('Beta_s_u_b / Beta_o_p_t');
        end
        xlim([min(data)-d max(data)+d]); ylim([0 1]); lineplot(0,'v','k--');
    end
end









