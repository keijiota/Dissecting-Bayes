clear all; close all;

load('dat_percep.mat');

xx = 0:0.01:1; result = [];
for condition = 1:2
    % fit subject average
    if condition == 1
        % 30, symmetric
        X = MeanCorrectProb(:,:,1)*2 ; Y = MeanEstProb(:,:,1) ;
        conname = '30 dots';
    else
        % 5, symmetric
        X = MeanCorrectProb(:,:,4)*2  ; Y = MeanEstProb(:,:,4) ;
        conname = '5 dots';
    end
    stdy = std(Y); x = X; y = Y;
    X = mean(X); Y = mean(Y); 
%     X = reshape(X, 1,153); Y = reshape(Y, 1,153); 
    
    %         model test
    %         a = 1; b = 0;
    %         Y = X.*a + b + (rand(1,9)*0.0001-0.00005);
    %         a = 0.8+rand*0.4; b = rand*0.2-0.1 ;
    %         Y = X.*a + b ;
    %         gamma = rand+0.5 ;
    %         Y = (X.^gamma ./ ((X.^gamma + (1-X).^gamma).^(1/gamma))); % Tversky & Kahneman 1992
    %         Y = exp(-(-log(X)).^gamma) ;
    
    % H0: P[S]_est = P[S] + e
    N = length(X) ; K = 0;
    a = 1 ; b = 0 ;
    Y_hut = X.*a + b ;
    Residual_h0(1,condition) = mean((Y-Y_hut).^2) ;
    nll_h0 = - length(X)/2 * (log(2*pi) + log(Residual_h0(1,condition)) + 1) ;% bigger is good
    AICc_h0(1,condition) = -2*(nll_h0) + 2*K + 2*(K*(K+1))/(N-K-1) ; % smaller is good
    y_h0 = xx.*a + b ;
    
    % H1: P[SU] = alpha * P[SL] + 0 + e
    N = length(X) ; K = 1;
    for rep = 1:5
        NLL = @(slope) [N/2 * (log(2*pi) + log(mean(((slope*X + 0) - Y).^2)) + 1)];
        [fitslope_pre mnll_pre] = fmincon(NLL,[rand*2],[],[],[],[],[],[]);
        if rep == 1
            fitslope = fitslope_pre; mnll = mnll_pre;
        elseif rep > 1 && mnll >= mnll_pre
            fitslope = fitslope_pre; mnll = mnll_pre;
        end
    end
    Y_hut = fitslope*X + 0;
    Slope(1,condition) = fitslope;
    Residual_h1(1,condition) = mean((Y-Y_hut).^2) ;
    nll_h1 = - length(X)/2 * (log(2*pi) + log(Residual_h1(1,condition)) + 1) ;
    AICc_h1(1,condition) = -2*(nll_h1) + 2*K + 2*(K*(K+1))/(N-K-1);
    y_h1 = xx*fitslope + 0;
    
    % H2: P[SU] = 1 * P[SL] + beta + e
    N = length(X) ; K = 1;
    for rep = 1:5
        NLL = @(intercept) [N/2 * (log(2*pi) + log(mean(((1*X + intercept) - Y).^2)) + 1)];
        [fitintercept_pre mnll_pre] = fmincon(NLL,[rand*2],[],[],[],[],[],[]);
        if rep == 1
            fitintercept = fitintercept_pre; mnll = mnll_pre;
        elseif rep > 1 && mnll >= mnll_pre
            fitintercept = fitintercept_pre; mnll = mnll_pre;
        end
    end
    Y_hut = 1*X + fitintercept;
    Intercept(1,condition) = fitintercept;
    Residual_h2(1,condition) = mean((Y-Y_hut).^2) ;
    nll_h2 = - length(X)/2 * (log(2*pi) + log(Residual_h2(1,condition)) + 1) ;
    AICc_h2(1,condition) = -2*(nll_h2) + 2*K + 2*(K*(K+1))/(N-K-1);
    y_h2 = 1*xx + fitintercept;
    
    % H1: P[S]_est = alpha * P[S] + beta + e
    N = length(X) ; K = 2;
    mdl = fitglm(X,Y) ;
    nll_h3 = mdl.LogLikelihood;
    Residual_h3(1,condition) = mean(mdl.Residuals.Raw.^2) ;
    alpha(1,condition) = mdl.Coefficients.Estimate(2,1) ;
    beta(1,condition) = mdl.Coefficients.Estimate(1,1) ;
    AICc_h3(1,condition) = -2*(nll_h3) + 2*K + 2*(K*(K+1))/(N-K-1);
    y_h3 = xx.*alpha(1,condition) + beta(1,condition) ;
    
    % H4: P[S]_est = P[S]^gamma / (P[S]^gamma + 1-P[S]^gamma)^(1/gamma), Tversky & Kahneman 1992
    N = length(X) ; K = 1;
    for rep = 1:5
        NLL = @(gamma) [N/2 * (log(2*pi) + log(mean(((X.^gamma./((X.^gamma + (1-X).^gamma).^(1/gamma))) - Y).^2)) + 1)] ;
        [fitgamma_pre mnll_pre] = fmincon(NLL,[rand*2],[],[],[],[],[],[]);
        if rep == 1
            fitgamma = fitgamma_pre; mnll = mnll_pre;
        elseif rep > 1 && mnll >= mnll_pre
            fitgamma = fitgamma_pre; mnll = mnll_pre;
        end
    end
    Y_hut = X.^fitgamma./((X.^fitgamma + (1-X).^fitgamma).^(1/fitgamma)) ;
    gamma_KT(1,condition) = fitgamma;
    Residual_h4(1,condition) = mean((Y-Y_hut).^2) ;
    nll_h4 = - length(X)/2 * (log(2*pi) + log(Residual_h4(1,condition)) + 1) ;
    AICc_h4(1,condition) = -2*(nll_h4) + 2*K + 2*(K*(K+1))/(N-K-1);
    y_h4 = xx.^fitgamma./((xx.^fitgamma + (1-xx).^fitgamma).^(1/fitgamma)) ;
    
    % H5: P[S]_est = exp(-(-ln(p)^gamma)),Prelec 1998
    N = length(X) ; K = 1;
    for rep = 1:5
        NLL = @(gamma) [N/2 * (log(2*pi) + log(mean(((exp(-(-log(X)).^gamma)) - Y).^2)) + 1)];
        [fitgamma_pre2 mnll_pre2] = fmincon(NLL,[rand*2],[],[],[],[],[],[]);
        if rep == 1
            fitgamma2 = fitgamma_pre2; mnll = mnll_pre2;
        elseif rep > 1 && mnll >= mnll_pre
            fitgamma2 = fitgamma_pre2; mnll = mnll_pre2;
        end
    end
    Y_hut = exp(-(-log(X)).^fitgamma2) ;
    gamma_Pr(1,condition) = fitgamma2;
    Residual_h5(1,condition) = mean((Y-Y_hut).^2) ;
    nll_h5 = - length(X)/2 * (log(2*pi) + log(Residual_h5(1,condition)) + 1) ;
    AICc_h5(1,condition) = -2*(nll_h5) + 2*K + 2*(K*(K+1))/(N-K-1);
    y_h5 = exp(-(-log(xx)).^fitgamma2) ;
    
    subplot(1,2,condition)
    plot(X,Y, 'ko', 'MarkerFaceColor','k') ; hold on
    errorbar(X, Y, stdy, 'ko', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', 0.1, 'CapSize',0) ;                
%     plot(xx, y_h0, 'r-');
%     plot(xx, y_h1, 'b-');
%     plot(xx, y_h2, 'g-');
%     plot(xx, y_h3, 'c-');
%     plot(xx, y_h4, 'y-');
%     plot(xx, y_h5, 'm-');
    plot(0:0.1:1,0:0.1:1,'k--');
    axis('square');
    xlim([0 1]); ylim([0 1]);
    xticks(0:0.2:1); yticks(0:0.2:1);
    xlabel('P[S]', 'FontName', 'Arial', 'FontSize', 10);
    ylabel('\pi(P[S])', 'FontName', 'Arial', 'FontSize', 10);
    title(strcat('Subject average, ',conname));
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);    
    
    for i = 1:9
        [h, p, ci, stat] = ttest(y(:,i), x(:,i),'alpha',0.05/9);
        meandiff = mean(y(:,i)) - mean(x(:,i));
        result = [result; condition, stat.tstat, stat.df, p, h, meandiff]; 
    end        
end

deltaAicc2null_30d = AICc_h0(:,1) - [AICc_h0(:,1) AICc_h1(:,1) AICc_h2(:,1) AICc_h3(:,1) AICc_h4(:,1) AICc_h5(:,1)];
% ngative means null is better
evidenceratio_30d = exp((deltaAicc2null_30d - deltaAicc2null_30d(1))/2) ;
% er of 0.33 means null is 3 times likely, er of 3 means null is 3 times less likely
deltaAicc2null_5d = AICc_h0(:,2) - [AICc_h0(:,2) AICc_h1(:,2) AICc_h2(:,2) AICc_h3(:,2) AICc_h4(:,2) AICc_h5(:,2)];
evidenceratio_5d = exp((deltaAicc2null_5d - deltaAicc2null_5d(1))/2) ;

result = [deltaAicc2null_30d([1,4:6])', deltaAicc2null_5d([1,4:6])', evidenceratio_30d([1,4:6])', evidenceratio_5d([1,4:6])']; 


C1 = [255 200 200]./255; C11 = [200 0 0]./255;
C2 = [179 179 255]./255; C22 = [25 25 255]./255;
C3 = [200 255 200]./255 ; C33 = [52 255 50]./255 ;
C4 = [204 247 255]./255 ; C44 = [80 180 255]./255 ;
C5 = [255 240 0]./255 ; C55 = [255 155 0]./255 ;
C6 = [240 210 240]./255 ; C66 = [255 0 255]./255 ;
ms = 8; cs = 10; ms2 = 10;

figure
subplot(1,2,1)
plot(1, AICc_h0(:,1), 's', 'MarkerEdgeColor', C11, 'MarkerFaceColor', C11, 'MarkerSize',ms2); hold on 
plot(2, AICc_h1(:,1), 's', 'MarkerEdgeColor', C22, 'MarkerFaceColor', C22, 'MarkerSize',ms2);
plot(3, AICc_h2(:,1), 's', 'MarkerEdgeColor', C33, 'MarkerFaceColor', C33, 'MarkerSize',ms2);
plot(4, AICc_h3(:,1), 's', 'MarkerEdgeColor', C44, 'MarkerFaceColor', C44, 'MarkerSize',ms2);
plot(5, AICc_h4(:,1), 's', 'MarkerEdgeColor', C55, 'MarkerFaceColor', C55, 'MarkerSize',ms2);
plot(6, AICc_h5(:,1), 's', 'MarkerEdgeColor', C66, 'MarkerFaceColor', C66, 'MarkerSize',ms2);
xlim([0.25 6.75]); xticks(1:1:6);
xticklabels({'H0(y~x)', 'H1(y~a*x+0)', 'H2(y~1*x+b)', 'H3(y~a*x+b)', 'H4(y~x^g/(x^g+(1-x)^g)^1^/^g', 'H5(y~exp(-(-ln(x)^g))'});
ylabel('AICc', 'FontName', 'Arial', 'FontSize', 10);
title('30 dots');
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 14, 'linewidth', 1.5,'TickLength',[0.025 0]);

subplot(1,2,2)
plot(1, AICc_h0(:,2), 's', 'MarkerEdgeColor', C11, 'MarkerFaceColor', C11, 'MarkerSize',ms2); hold on
plot(2, AICc_h1(:,2), 's', 'MarkerEdgeColor', C22, 'MarkerFaceColor', C22, 'MarkerSize',ms2);
plot(3, AICc_h2(:,2), 's', 'MarkerEdgeColor', C33, 'MarkerFaceColor', C33, 'MarkerSize',ms2);
plot(4, AICc_h3(:,2), 's', 'MarkerEdgeColor', C44, 'MarkerFaceColor', C44, 'MarkerSize',ms2);
plot(5, AICc_h4(:,2), 's', 'MarkerEdgeColor', C55, 'MarkerFaceColor', C55, 'MarkerSize',ms2);
plot(6, AICc_h5(:,2), 's', 'MarkerEdgeColor', C66, 'MarkerFaceColor', C66, 'MarkerSize',ms2);
xlim([0.25 6.75]); xticks(1:1:6);
xticklabels({'H0(y~x)', 'H1(y~a*x+0)', 'H2(y~1*x+b)', 'H3(y~a*x+b)', 'H4(y~x^g/(x^g+(1-x)^g)^1^/^g', 'H5(y~exp(-(-ln(x)^g))'});
ylabel('AICc', 'FontName', 'Arial', 'FontSize', 10);
title('5 dots');
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 14, 'linewidth', 1.5,'TickLength',[0.025 0]);

