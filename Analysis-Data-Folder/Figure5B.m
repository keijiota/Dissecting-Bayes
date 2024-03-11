clear; close all;

load('dat_percep.mat');

xx = 0:0.01:1; result = [];
colormap = [0 0 207; 0 111 255; 0  255 255; 111 255 143;
    207 255 48; 255 191 0; 255 96 0; 200 0  0; 100  0  0];
for condition = 1:2
    if condition == 1
        % 30, symmetric
        x = MeanEstProb(:,:,1) ; y1 = MeanEstProb(:,:,2); y2 = MeanEstProb(:,:,3);
        conname = '30-point';
    else
        % 5, symmetric
        x = MeanEstProb(:,:,4) ; y1 = MeanEstProb(:,:,5); y2 = MeanEstProb(:,:,6);
        conname = '5-point';
    end
    X = mean(x); Y = mean(y1+y2);
    
    % H0: P[S]_est = P[S] + e
    N = length(X) ; K = 0;
    a = 1 ; b = 0 ;
    Y_hut = X.*a + b ;
    Residual_h0(1,condition) = mean((Y-Y_hut).^2) ;
    nll_h0 = - length(X)/2 * (log(2*pi) + log(Residual_h0(1,condition)) + 1) ;% bigger is good
    AICc_h0(1,condition) = -2*(nll_h0) + 2*K + 2*(K*(K+1))/(N-K-1) ; % smaller is good
    BIC_h0(1,condition) = -2*(nll_h0) + log(N)*K;
    
    y_h0 = xx.*a + b ;
    
    % H1: P[SU] = P[SL] + b + e (b > 0)
    N = length(X) ; K = 1;
    for rep = 1:5
        NLL = @(params) [N/2 * (log(2*pi) + log(mean(((params(1)*X + params(2)) - Y).^2)) + 1)];
        [fitparams_pre mnll_pre] = fmincon(NLL,[rand(2,1)*2],[],[],[],[],[1 0],[1 inf]);
        if rep == 1
            fitparams = fitparams_pre; mnll = mnll_pre;
        elseif rep > 1 && mnll >= mnll_pre
            fitparams = fitparams_pre; mnll = mnll_pre;
        end
    end
    Y_hut = fitparams(1)*X + fitparams(2);
    Alpha(1,condition) = fitparams(1);
    Beta(1,condition) = fitparams(2);
    Residual_h1(1,condition) = mean((Y-Y_hut).^2) ;
    nll_h1 = - length(X)/2 * (log(2*pi) + log(Residual_h1(1,condition)) + 1) ;
    AICc_h1(1,condition) = -2*(nll_h1) + 2*K + 2*(K*(K+1))/(N-K-1);
    BIC_h1(1,condition) = -2*(nll_h1) + log(N)*K;    
    y_h1 = xx*fitparams(1) + fitparams(2);
    
    b(1,condition) = fitparams(2);
    % H2: P[SU] = P[SL] + b + e (b < 0)
    clear fitparams fitparams_pre mnll mnll_pre
    N = length(X) ; K = 1;
    for rep = 1:5
        NLL = @(params) [N/2 * (log(2*pi) + log(mean(((params(1)*X + params(2)) - Y).^2)) + 1)];
        [fitparams_pre mnll_pre] = fmincon(NLL,[rand(2,1)*2],[],[],[],[],[1 -inf],[1 0]);
        if rep == 1
            fitparams = fitparams_pre; mnll = mnll_pre;
        elseif rep > 1 && mnll >= mnll_pre
            fitparams = fitparams_pre; mnll = mnll_pre;
        end
    end
    Y_hut = fitparams(1)*X + fitparams(2);
    Alpha2(1,condition) = fitparams(1);
    Beta2(1,condition) = fitparams(2);
    Residual_h2(1,condition) = mean((Y-Y_hut).^2) ;
    nll_h2 = - length(X)/2 * (log(2*pi) + log(Residual_h2(1,condition)) + 1) ;
    AICc_h2(1,condition) = -2*(nll_h2) + 2*K + 2*(K*(K+1))/(N-K-1);
    BIC_h2(1,condition) = -2*(nll_h2) + log(N)*K;
    y_h2 = xx*fitparams(1) + fitparams(2);

    subplot(1,2,condition)
    for i = 1:9
        errorbar(mean(x(:,i)), mean([y1(:,i)+y2(:,i)]), 2*std(x(:,i))/sqrt(17), 'horizontal', 'k', 'linewidth', 0.5, 'MarkerSize', 0.1, 'CapSize',5); hold on
        errorbar(mean(x(:,i)), mean([y1(:,i)+y2(:,i)]), [2*std(y1(:,i))/sqrt(17) + 2*std(y2(:,i))/sqrt(17)]/2, 'k', 'linewidth', 0.5, 'MarkerSize', 0.1, 'CapSize',5); hold on
        plot(mean(x(:,i)), mean([y1(:,i)+y2(:,i)]), 'ko', 'MarkerFaceColor',[colormap(i,:)/255], 'Markersize',7) ; hold on
    end

    plot(0:0.1:1.2,0:0.1:1.2,'k--');
    plot(xx, y_h1,'k-');
    lineplot(1,'h', 'k:'); lineplot(1,'v', 'k:');
    axis('square');
    xlim([0 1.05]); ylim([0 1.05]);
    xticks(0:0.2:1.2); yticks(0:0.2:1.2);
    xlabel('\pi(P[S])', 'FontName', 'Arial', 'FontSize', 10);
    ylabel('\pi(P[SU])+\pi(P[SL])', 'FontName', 'Arial', 'FontSize', 10);
    title(conname);
    set(gca, 'Fontname', 'Arial', 'Fontsize', 14, 'linewidth', 1.5, 'TickLength',[0.025 0]);
    
    for i = 1:9
        [h, p, ci, stat] = ttest(y1(:,i)+y2(:,i), x(:,i),'alpha',0.05/9);
        meandiff = mean(y1(:,i)+y2(:,i)) - mean(x(:,i));
        result = [result; condition, stat.tstat, stat.df, p, h, meandiff];
    end
    
    Meandiff(:,condition) = mean(y1+y2) - mean(x);
    
end

deltaAicc2null_30d = AICc_h0(:,1) - [AICc_h0(:,1) AICc_h1(:,1) AICc_h2(:,1)];
% negative means null is better
evidenceratio_30d = exp((deltaAicc2null_30d - deltaAicc2null_30d(1))/2) ;
deltaBIC2null_30d = BIC_h0(:,1) - [BIC_h0(:,1) BIC_h1(:,1) BIC_h2(:,1)];


% er of 0.33 means null is 3 times likely, er of 3 means null is 3 times less likely
deltaAicc2null_5d = AICc_h0(:,2) - [AICc_h0(:,2) AICc_h1(:,2) AICc_h2(:,2)] ;
evidenceratio_5d = exp((deltaAicc2null_5d - deltaAicc2null_5d(1))/2) ;
deltaBIC2null_5d = BIC_h0(:,2) - [BIC_h0(:,2) BIC_h1(:,2) BIC_h2(:,2)];

result = [deltaBIC2null_30d' deltaBIC2null_5d', deltaAicc2null_30d', deltaAicc2null_5d', evidenceratio_30d', evidenceratio_5d']; 

beta = [Beta; Beta2]
