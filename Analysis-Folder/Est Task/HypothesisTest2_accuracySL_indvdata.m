clear all; close all;

load('dat_percep.mat');

xx = 0:0.01:1;
for condition = 1:2
    for subi = 1:17 % fit individual data
        if condition == 1
            % 30, symmetric
            X = MeanCorrectProb(subi,:,1) ; 
            Y = MeanEstProb(subi,:,3) ; % SL            
            figure(1)
            conname = '30 dots';
        else
            % 5, symmetric
            X = MeanCorrectProb(subi,:,4)  ; 
            Y = MeanEstProb(subi,:,6) ; % SL            
            figure(2)
            conname = '5 dots';
        end        
        
        Y(Y==0) = 0.00000001;
        
        % H0: linear log odds transformation with slope of 1
        N = length(X) ; K = 0;
        a = 1 ; b = 0 ;
        Y_hut = X.*a + b ;
        Residual_h0(subi,condition) = mean((Y-Y_hut).^2) ;
        nll_h0 = - length(X)/2 * (log(2*pi) + log(Residual_h0(subi,condition)) + 1) ;% bigger is good
        AICc_h0(subi,condition) = -2*(nll_h0) + 2*K + 2*(K*(K+1))/(N-K-1) ; % smaller is good
        y_h0 = xx.*a + b ;

        % H1: linear log odds transformation
        N = length(X) ; K = 2;
        logodds = @(p) log(p./(1-p));
        Lo_Pi_p = logodds(Y);
        Lo_p = logodds(X);
        
        fun = @(pm)  -sum(-log(pm(3)) - (Lo_Pi_p - (pm(1) * logodds(X) +(1-pm(1)) * logodds(pm(2)))).^2/(2*(pm(3).^2)));
        pm = fmincon(fun,[rand(1,3)],[],[],[],[],[0 0 0],[100 1 100]) ;                
        
        LoPip_hut = pm(1) * logodds(X) + (1-pm(1)) * logodds(pm(2));
        Y_hut = exp(LoPip_hut) ./ (1 + exp(LoPip_hut)) ; % inverse of log odds
        Gamma(subi,condition) = pm(1);
        CrossOver(subi,condition) = pm(2);
        Residual_h1(subi,condition) = mean((Y - Y_hut).^2) ;
        nll_h1 = - length(X)/2 * (log(2*pi) + log(Residual_h1(subi,condition)) + 1) ;% bigger is good
        AICc_h1(subi,condition) = -2*(nll_h1) + 2*K + 2*(K*(K+1))/(N-K-1) ; % smaller is good
        y_h1 = pm(1) * logodds(xx) + (1-pm(1)) * logodds(pm(2));
        y_h1 = exp(y_h1) ./ (1 + exp(y_h1)) ;
                
        % H2: P[S]_est = P[S]^gamma / (P[S]^gamma + 1-P[S]^gamma)^(1/gamma), Tversky & Kahneman 1992
        N = length(X) ; K = 1;
        fun = @(pm)  -sum(-log(pm(2)) - (Y - (X.^pm(1)./((X.^pm(1) + (1-X).^pm(1)).^(1/pm(1))))).^2/(2*(pm(2).^2)));
        pm = fmincon(fun, [rand(1,2)],[],[],[],[],[0 0],[100 100]);
        Gamma_kt(subi,condition) = pm(1);
        Y_hut = X.^pm(1)./((X.^pm(1) + (1-X).^pm(1)).^(1/pm(1)));
        Residual_h2(subi,condition) = mean((Y-Y_hut).^2);
        nll_h2 = - length(X)/2 * (log(2*pi) + log(Residual_h2(subi,condition)) + 1) ;
        AICc_h2(subi,condition) = -2*(nll_h2) + 2*K + 2*(K*(K+1))/(N-K-1) ;
        
        % H3: P[S]_est = exp(-(-ln(p)^gamma)),Prelec 1998
        N = length(X) ; K = 1;
        fun = @(pm)  -sum(-log(pm(2)) - (Y - (exp(-(-log(X)).^pm(1)))).^2/(2*(pm(2).^2)));
        pm = fmincon(fun, [rand(1,2)],[],[],[],[],[0 0],[100 100]);
        Gamma_pr(subi,condition) = pm(1);
        Y_hut = exp(-(-log(X)).^pm(1));
        Residual_h3(subi,condition) = mean((Y-Y_hut).^2) ;
        nll_h3 = - length(X)/2 * (log(2*pi) + log(Residual_h3(subi,condition)) + 1) ;
        AICc_h3(subi,condition) = -2*(nll_h3) + 2*K + 2*(K*(K+1))/(N-K-1) ;
        
        
        subplot(3,6,subi)
        plot(X,Y, 'ko', 'MarkerFaceColor','k','MarkerSize',9) ; hold on
        plot(0:0.1:1,0:0.1:1,'k--'); 
        plot(xx, y_h1, 'k-');
        axis('square');        
        xlim([0 1]); ylim([0 1]);
        xticks(0:0.5:1); yticks(0:0.5:1);
        xlabel('P[S]', 'FontName', 'Arial', 'FontSize', 10);
        ylabel('\pi(P[S])', 'FontName', 'Arial', 'FontSize', 10);
        title(strcat('Sub', num2str(subi),',',conname));
        set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);                
    end
end

deltaAicc2null_30d = AICc_h0(:,1) - [AICc_h0(:,1) AICc_h1(:,1)]; % negative means null is better
evidenceratio_30d = exp((deltaAicc2null_30d - deltaAicc2null_30d(1))/2) ;
% er of 0.33 means null is 3 times likely, er of 3 means null is 3 times less likely
deltaAicc2null_5d = AICc_h0(:,2) - [AICc_h0(:,2) AICc_h1(:,2)]; 
evidenceratio_5d = exp((deltaAicc2null_5d - deltaAicc2null_5d(1))/2) ;

meandeltaAicc2null_30d = mean(AICc_h0(:,1)) - [mean(AICc_h0(:,1)) mean(AICc_h1(:,1))]; 
meanevidenceratio_30d = exp((meandeltaAicc2null_30d - meandeltaAicc2null_30d(1))/2) ;
meandeltaAicc2null_5d = mean(AICc_h0(:,2)) - [mean(AICc_h0(:,2)) mean(AICc_h1(:,2))];
meanevidenceratio_5d = exp((meandeltaAicc2null_5d - meandeltaAicc2null_5d(1))/2) ;

result = [meandeltaAicc2null_30d([1:2])', meandeltaAicc2null_5d([1:2])', meanevidenceratio_30d([1:2])', meanevidenceratio_5d([1:2])'] 



