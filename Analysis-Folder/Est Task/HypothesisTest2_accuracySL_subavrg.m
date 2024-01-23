clear all; close all;

load('dat_percep.mat');

% analysis of SU and SL 

xx = 0:0.01:1; result = [];
for condition = 1:2
    % fit subject average
    if condition == 1
        % 30, symmetric
        X = MeanCorrectProb(:,:,1);   
        Y = MeanEstProb(:,:,3) ; % SL        
        conname = '30 dots';
    else
        % 5, symmetric
        X = MeanCorrectProb(:,:,4); 
        Y = MeanEstProb(:,:,6); % SL        
        conname = '5 dots';
    end
    stdy = std(Y); x = X; y = Y;
    X = mean(X); Y = mean(Y); 
    
    % H0: slope of 1
    N = length(X) ; K = 0;
    a = 1 ; b = 0 ;
    Y_hut = X.*a + b ;
    Residual_h0(1,condition) = mean((Y-Y_hut).^2) ;
    nll_h0 = - length(X)/2 * (log(2*pi) + log(Residual_h0(1,condition)) + 1) ;% bigger is good
    AICc_h0(1,condition) = -2*(nll_h0) + 2*K + 2*(K*(K+1))/(N-K-1) ; % smaller is good
    BIC_h0(1,condition) = -2*(nll_h0) + log(N)*K;
    y_h0 = xx.*a + b ;        
    
    % H1: linear log odds transformation
    N = length(X) ; K = 2;         
    logodds = @(p) log(p./(1-p));          
    Lo_Pi_p = logodds(Y); 
    Lo_p = logodds(X); 
    
    fun = @(prm)  -sum(-log(prm(3)) - (Lo_Pi_p - (prm(1) * logodds(X) +(1-prm(1)) * logodds(prm(2)))).^2/(2*(prm(3).^2)));    
    pm = fmincon(fun,[rand(1,3)]);
    
    LoPip_hut = pm(1) * logodds(X) + (1-pm(1)) * logodds(pm(2)); 
    Y_hut = exp(LoPip_hut) ./ (1 + exp(LoPip_hut)) ; % inverse of log odds
    Gamma(1,condition) = pm(1);
    CrossOver(1,condition) = pm(2);
    Residual_h1(1,condition) = mean((Y - Y_hut).^2) ;
    nll_h1 = - length(X)/2 * (log(2*pi) + log(Residual_h1(1,condition)) + 1) ;% bigger is good
    AICc_h1(1,condition) = -2*(nll_h1) + 2*K + 2*(K*(K+1))/(N-K-1) ; % smaller is good
    BIC_h1(1,condition) = -2*(nll_h1) + log(N)*K;

    y_h1 = pm(1) * logodds(xx) + (1-pm(1)) * logodds(pm(2));    
    y_h1 = exp(y_h1) ./ (1 + exp(y_h1)) ; 
        
    % H2: P[S]_est = P[S]^gamma / (P[S]^gamma + 1-P[S]^gamma)^(1/gamma), Tversky & Kahneman 1992
    N = length(X) ; K = 1;    
    fun = @(pm)  -sum(-log(pm(2)) - (Y - (X.^pm(1)./((X.^pm(1) + (1-X).^pm(1)).^(1/pm(1))))).^2/(2*(pm(2).^2)));    
    pm = fmincon(fun, [rand(1,2)]);    
    Gamma_kt(1,condition) = pm(1);    
    Y_hut = X.^pm(1)./((X.^pm(1) + (1-X).^pm(1)).^(1/pm(1)));     
    Residual_h2(1,condition) = mean((Y-Y_hut).^2);    
    nll_h2 = - length(X)/2 * (log(2*pi) + log(Residual_h2(1,condition)) + 1) ; 
    AICc_h2(1,condition) = -2*(nll_h2) + 2*K + 2*(K*(K+1))/(N-K-1) ; 
    BIC_h2(1,condition) = -2*(nll_h2) + log(N)*K;

    % H3: P[S]_est = exp(-(-ln(p)^gamma)),Prelec 1998
    N = length(X) ; K = 1;    
    fun = @(pm)  -sum(-log(pm(2)) - (Y - (exp(-(-log(X)).^pm(1)))).^2/(2*(pm(2).^2)));    
    pm = fmincon(fun, [rand(1,2)]);    
    Gamma_pr(1,condition) = pm(1);    
    Y_hut = exp(-(-log(X)).^pm(1));
    Residual_h3(1,condition) = mean((Y-Y_hut).^2) ;    
    nll_h3 = - length(X)/2 * (log(2*pi) + log(Residual_h3(1,condition)) + 1) ; 
    AICc_h3(1,condition) = -2*(nll_h3) + 2*K + 2*(K*(K+1))/(N-K-1) ;         
    BIC_h3(1,condition) = -2*(nll_h3) + log(N)*K;

    figure(1)
    subplot(1,2,condition)
    plot(X,y, 'ko', 'MarkerFaceColor','w','MarkerSize', 3) ; hold on    
    plot(X,Y, 'ko', 'MarkerFaceColor','k') ; hold on
%     errorbar(X, Y, stdy, 'ko', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', 0.1, 'CapSize',0) ;                
    plot(0:0.1:1,0:0.1:1,'k--');
    plot(xx, y_h1, 'k-');
    axis('square');
    xlim([0 0.75]); ylim([0 0.75]);
    xticks(0:0.1:1); yticks(0:0.1:1);
    xlabel('P[SL]', 'FontName', 'Arial', 'FontSize', 10);
    ylabel('\pi(P[SL])', 'FontName', 'Arial', 'FontSize', 10);
    title(strcat('Subject average, ',conname));
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);    
end


deltaAicc2null_30d = AICc_h0(:,1) - [AICc_h0(:,1) AICc_h1(:,1)]; % ngative means null is better
evidenceratio_30d = exp((deltaAicc2null_30d - deltaAicc2null_30d(1))/2) ;

% er of 0.33 means null is 3 times likely, er of 3 means null is 3 times less likely
deltaAicc2null_5d = AICc_h0(:,2) - [AICc_h0(:,2) AICc_h1(:,2)];
evidenceratio_5d = exp((deltaAicc2null_5d - deltaAicc2null_5d(1))/2) ;

result = [deltaAicc2null_30d([1:2])', deltaAicc2null_5d([1:2])', evidenceratio_30d([1:2])', evidenceratio_5d([1:2])'] 

deltaAicc2h1_30d =  AICc_h1(:,1) - [AICc_h1(:,1) AICc_h2(:,1) AICc_h3(:,1)] ;
deltaAicc2h1_5d =  AICc_h1(:,2) - [AICc_h1(:,2) AICc_h2(:,2) AICc_h3(:,2)] ;
evidenceratioh1_30d = exp((deltaAicc2h1_30d - deltaAicc2h1_30d(1))/2) ;
evidenceratioh1_5d = exp((deltaAicc2h1_5d - deltaAicc2h1_5d(1))/2) ;

result2 = [deltaAicc2h1_30d', deltaAicc2h1_5d', evidenceratioh1_30d', evidenceratioh1_5d'] 

CrossOver
gamma = [Gamma; Gamma_kt; Gamma_pr]


deltaAicc2h1_30d =  AICc_h0(:,1) - [AICc_h0(:,1) AICc_h1(:,1) AICc_h2(:,1) AICc_h3(:,1)] ;
deltaAicc2h1_5d =  AICc_h0(:,2) - [AICc_h0(:,2) AICc_h1(:,2) AICc_h2(:,2) AICc_h3(:,2)] ;
evidenceratioh1_30d = exp((deltaAicc2h1_30d - deltaAicc2h1_30d(1))/2) ;
evidenceratioh1_5d = exp((deltaAicc2h1_5d - deltaAicc2h1_5d(1))/2) ;

deltaBIC2h1_30d =  BIC_h0(:,1) - [BIC_h0(:,1) BIC_h1(:,1) BIC_h2(:,1) BIC_h3(:,1)] ;
deltaBIC2h1_5d =  BIC_h0(:,2) - [BIC_h0(:,2) BIC_h1(:,2) BIC_h2(:,2) BIC_h3(:,2)] ;

result3 = [deltaBIC2h1_30d' deltaBIC2h1_5d' deltaAicc2h1_30d', deltaAicc2h1_5d', evidenceratioh1_30d', evidenceratioh1_5d'] 
