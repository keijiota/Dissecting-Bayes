clear; close all ;

load regression_partial;
load sim_weight;

xt = {'Pmin','P50%','Pmax'};

C1 = [150 150 150]./255 ; C11 = [10 10 10]./255 ;
C2 = [255 200 200]./255 ; C22 = [255 75 75]./255 ;
ms = 8; ms2 = 5; cs = 10;

% ms = 7; xmin = [-1.5 -1.5 -1.5]; xmax = [1.5 1.5 1.5];
tn = {'30 dots, 0 penalty', '5 dots, 0 penalty', '30 dots, -500 penalty','5 dots, -500 penalty'};
nsub = 17; Con = 1:4;
b_opt = meanW;
b_opt{1,1} = b_opt{1,1}([1,15,30]); b_opt{1,2} = b_opt{1,2}([1,15,30]);
b_opt{1,3} = b_opt{1,3}([1,3,5]);   b_opt{1,4} = b_opt{1,4}([1,3,5]);

%% plot beta
figure; ic = 0;
for con = [1 3 2 4]
    data = coefficient_obs(:,:,con);
    regopt = coefficient_opt(:,:,con);
    nfeatures = size(coefficient_obs(:,:,con),2); ic = ic + 1;
    
    subplot(2,2,Con(ic)); hold on
    avbetaobs = mean(data,1); sebetaobs = std(data,[],1); % / sqrt(nsub) ;
    p1 = plot(1:nfeatures, b_opt{1,con}, '-', 'Color','k','linewidth', 1);
    p2 = plot(1:nfeatures, mean(regopt,1), '--s', 'Color','k','MarkerFaceColor','k','linewidth', 1);    
    p3 = plot(1:nfeatures, avbetaobs, 's:', 'color', C22, 'MarkerEdgeColor', C22,'MarkerFaceColor', C22,'linewidth', 1.5, 'MarkerSize', ms2);
%     p4 = plot(1:nfeatures, data, '-', 'Color',C2);
    p4 = plot(1:nfeatures, data, 'wo', 'MarkerFaceColor',C2, 'MarkerSize',5);
    p5 = errorbar(1:nfeatures, avbetaobs, sebetaobs, 's:', 'color', C22, 'MarkerFaceColor', C22,'linewidth', 1.5, 'MarkerSize', ms2, 'CapSize',cs);
    
    ylim([-1 0.75]); yticks(-1.5:0.5:1.5); title(tn(ic)); myfig;
    ylabel('Beta'); xlim([0.75 nfeatures+0.25]); xticklabels(xt);
end
legend([p1 p2 p3], {'Simulation','Normal regression to S_i_d_e_a_l','Normal regression to S_s_u_b'});


%% plot ratio
for i = 1:2
figure; ic = 0;
for con = [1 3 2 4]
    nfeatures = size(coefficient_obs(:,:,con),2); ic = ic + 1;
    
    if i == 1
        data = coefficient_obs(:,:,con) ./ b_opt{1,con};
    else
        data = coefficient_obs(:,:,con) ./ coefficient_opt(:,:,con);
    end
    Betaratio{i,con} = data ;
    
    subplot(2,2,Con(ic)); hold on
    ylabel('Beta_s_u_b / Beta_o_p_t');
    title(tn(con));  myfig ;
    ylim([-5 10]); yticks(-100:5:100);
    xlim([0.75 nfeatures+0.25]); xticklabels(xt);
    lineplot(0, 'h','k--'); lineplot(1,'h', 'k--');    
    
    avbetaratio = mean(data,1); sebetaratio = std(data,[],1); % / sqrt(nsub) ;
    p3 = plot(1:nfeatures, data, 'wo', 'MarkerFaceColor',C2, 'MarkerSize',5);    
    p4 = errorbar(1:nfeatures, avbetaratio, sebetaratio, 's:', 'color', C22, 'MarkerFaceColor', C22,'linewidth', 1.5, 'MarkerSize', ms2, 'CapSize',cs);
    
end
    legend(p4, {'Normal regression to S_s_u_b'});
end




