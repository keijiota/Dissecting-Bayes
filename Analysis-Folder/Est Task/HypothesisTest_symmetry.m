clear all; close all;

load('dat_percep.mat');

for condition = 1:2
    for subi = 1:17
        if condition == 1
            % 30, asymmetric(upper) and asymmetric(lower)
            X = MeanEstProb(subi,:,2) ; Y = MeanEstProb(subi,:,3) ;
            figure(1)
            conname = '30 dots';
        else
            % 5, asymmetric(upper) and asymmetric(lower)
            X = MeanEstProb(subi,:,5) ; Y = MeanEstProb(subi,:,6) ;
            figure(2)
            conname = '5 dots';
        end        
        xx = 0:0.01:1;
        
%         model test
%         a = 1; b = 0;
%         Y = X.*a + b + (rand(1,9)*0.1-0.05);
%         a = 0.8+rand*0.4; b = rand*0.2-0.1 ;
%         Y = X.*a + b ; 
%         Y = X.*0.9 + 0 + (rand(1,9)*0.01-0.005);

        % H0: P[SU] = P[SL] + e
        N = length(X) ; K = 0;
        a = 1 ; b = 0 ;
        Y_hut = X.*a + b ;
        Residual_h0(subi,condition) = mean((Y-Y_hut).^2) ;
        nll_h0 = - length(X)/2 * (log(2*pi) + log(Residual_h0(subi,condition)) + 1) ;% bigger is good
        AICc_h0(subi,condition) = -2*(nll_h0) + 2*K + 2*(K*(K+1))/(N-K-1) ; % smaller is good
        y_h0 = xx.*a + b ; 
       
        % H1: P[SU] = alpha * P[SL] + 0 + e
        N = length(X) ; K = 1;
         for rep = 1:5
          NLL = @(slope) [N/2 * (log(2*pi) + log(mean(( (slope*X + 0) - Y).^2)) + 1)];
          [fitslope_pre mnll_pre] = fmincon(NLL,[rand*2],[],[],[],[],[],[]);
            if rep == 1
                fitslope = fitslope_pre; mnll = mnll_pre;
            elseif rep > 1 && mnll >= mnll_pre
                fitslope = fitslope_pre; mnll = mnll_pre;                
            end
         end
        Y_hut = fitslope*X + 0;
        Slope(subi,condition) = fitslope; 
        Residual_h1(subi,condition) = mean((Y-Y_hut).^2) ;
        nll_h1 = - length(X)/2 * (log(2*pi) + log(Residual_h1(subi,condition)) + 1) ;
        AICc_h1(subi,condition) = -2*(nll_h1) + 2*K + 2*(K*(K+1))/(N-K-1);
        y_h1 = xx*fitslope + 0; 
        
        % H2: P[SU] = 1 * P[SL] + beta + e
        N = length(X) ; K = 1;
         for rep = 1:5
          NLL = @(intercept) [N/2 * (log(2*pi) + log(mean(( (1*X + intercept) - Y).^2)) + 1)];
          [fitintercept_pre mnll_pre] = fmincon(NLL,[rand*2],[],[],[],[],[],[]);
            if rep == 1
                fitintercept = fitintercept_pre; mnll = mnll_pre;
            elseif rep > 1 && mnll >= mnll_pre
                fitintercept = fitintercept_pre; mnll = mnll_pre;                
            end
         end
        Y_hut = 1*X + fitintercept;
        Intercept(subi,condition) = fitintercept; 
        Residual_h2(subi,condition) = mean((Y-Y_hut).^2) ;
        nll_h2 = - length(X)/2 * (log(2*pi) + log(Residual_h2(subi,condition)) + 1) ;
        AICc_h2(subi,condition) = -2*(nll_h2) + 2*K + 2*(K*(K+1))/(N-K-1);
        y_h2 = xx*fitslope + 0;         
       
        % H3: P[SU] = alpha * P[SL] + beta + e
        N = length(X) ; K = 2;
        mdl = fitglm(X,Y) ;
        nll_h1 = mdl.LogLikelihood;
        Residual_h3(subi,condition) = mean(mdl.Residuals.Raw.^2) ;
        alpha(subi,condition) = mdl.Coefficients.Estimate(2,1) ;
        beta(subi,condition) = mdl.Coefficients.Estimate(1,1) ;
        AICc_h3(subi,condition) = -2*(nll_h1) + 2*K + 2*(K*(K+1))/(N-K-1);
        y_h3 = xx.*alpha(subi,condition) + beta(subi,condition) ;
        
        subplot(3,6,subi)
        plot(X,Y, 'ko', 'MarkerFaceColor','k') ; hold on
        plot(xx, y_h0, 'r-', 'linewidth',1.5);
        plot(xx, y_h1, 'b-', 'linewidth',1.5);
        plot(xx, y_h2, 'g-', 'linewidth',1.5);
        plot(xx, y_h3, 'c-', 'linewidth',1.5);        
        axis('square');
        xlim([0 0.7]); ylim([0 0.7]); 
        xticks(0:0.1:0.7); yticks(0:0.1:0.7);
        xlabel('\pi(P[SU])', 'FontName', 'TimesNewRomans', 'FontSize', 10);
        ylabel('\pi(P[SL])', 'FontName', 'Arial', 'FontSize', 10);
        title(strcat('Sub ', num2str(subi),', ',conname));
        set(gca, 'Fontname', 'Arial', 'Fontsize', 14, 'linewidth', 1.5, 'TickLength',[0.025 0]);

    end
end

deltaAicc2null_30d = [AICc_h0(:,1) AICc_h1(:,1) AICc_h2(:,1) AICc_h3(:,1)] - AICc_h0(:,1);
% positive means null is better
evidenceratio_30d = exp((deltaAicc2null_30d - deltaAicc2null_30d(1))/2) ;
% er of 3 means null is 3 times likely, er of 0.33 means null is 3 times less likely
deltaAicc2null_5d = [AICc_h0(:,2) AICc_h1(:,2) AICc_h2(:,2) AICc_h3(:,2)] - AICc_h0(:,2);
evidenceratio_5d = exp((deltaAicc2null_5d - deltaAicc2null_5d(1))/2) ;

meandeltaAicc2null_30d = [mean(AICc_h0(:,1)) mean(AICc_h1(:,1)) mean(AICc_h2(:,1)) mean(AICc_h3(:,1))] - mean(AICc_h0(:,1));
meanevidenceratio_30d = exp((meandeltaAicc2null_30d - meandeltaAicc2null_30d(1))/2) ; 
meandeltaAicc2null_5d = [mean(AICc_h0(:,2)) mean(AICc_h1(:,2)) mean(AICc_h2(:,2)) mean(AICc_h3(:,2))] - mean(AICc_h0(:,2));
meanevidenceratio_5d = exp((meandeltaAicc2null_5d - meandeltaAicc2null_5d(1))/2) ; 


C1 = [255 200 200]./255; C11 = [200 0 0]./255; 
C2 = [179 179 255]./255; C22 = [25 25 255]./255;
C3 = [200 255 200]./255 ; C33 = [52 255 50]./255 ;
C4 = [204 247 255]./255 ; C44 = [80 180 255]./255 ; 
ms = 8; cs = 10; ms2 = 10;

figure
subplot(1,2,1)
plot(1, AICc_h0(:,1), 'o', 'MarkerEdgeColor', C1, 'MarkerFaceColor', C1, 'MarkerSize',ms); hold on
plot(2, AICc_h1(:,1), 'o', 'MarkerEdgeColor', C2, 'MarkerFaceColor', C2, 'MarkerSize',ms);
plot(3, AICc_h2(:,1), 'o', 'MarkerEdgeColor', C3, 'MarkerFaceColor', C3, 'MarkerSize',ms);
plot(4, AICc_h3(:,1), 'o', 'MarkerEdgeColor', C4, 'MarkerFaceColor', C4, 'MarkerSize',ms);
plot(1, mean(AICc_h0(:,1)), 's', 'MarkerEdgeColor', C11, 'MarkerFaceColor', C11, 'MarkerSize',ms2);
plot(2, mean(AICc_h1(:,1)), 's', 'MarkerEdgeColor', C22, 'MarkerFaceColor', C22, 'MarkerSize',ms2);
plot(3, mean(AICc_h2(:,1)), 's', 'MarkerEdgeColor', C33, 'MarkerFaceColor', C33, 'MarkerSize',ms2);
plot(4, mean(AICc_h3(:,1)), 's', 'MarkerEdgeColor', C44, 'MarkerFaceColor', C44, 'MarkerSize',ms2);
errorbar(1, mean(AICc_h0(:,1)), std(AICc_h0(:,1)), 'o', 'color', C11, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
errorbar(2, mean(AICc_h1(:,1)), std(AICc_h1(:,1)), 'o', 'color', C22, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
errorbar(3, mean(AICc_h2(:,1)), std(AICc_h2(:,1)), 'o', 'color', C33, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
errorbar(4, mean(AICc_h3(:,1)), std(AICc_h3(:,1)), 'o', 'color', C44, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
xlim([0.25 4.75]); xticks(1:1:4);
xticklabels({'H0(y~x)', 'H1(y~a*x+0)', 'H2(y~1*x+b)', 'H3(y~a*x+b)'});
ylabel('AICc', 'FontName', 'Arial', 'FontSize', 10);
title('30 dots');
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 14, 'linewidth', 1.5,'TickLength',[0.025 0]);

subplot(1,2,2)
plot(1, AICc_h0(:,2), 'o', 'MarkerEdgeColor', C1, 'MarkerFaceColor', C1, 'MarkerSize',ms); hold on
plot(2, AICc_h1(:,2), 'o', 'MarkerEdgeColor', C2, 'MarkerFaceColor', C2, 'MarkerSize',ms);
plot(3, AICc_h2(:,2), 'o', 'MarkerEdgeColor', C3, 'MarkerFaceColor', C3, 'MarkerSize',ms); hold on
plot(4, AICc_h3(:,2), 'o', 'MarkerEdgeColor', C4, 'MarkerFaceColor', C4, 'MarkerSize',ms);
plot(1, mean(AICc_h0(:,2)), 's', 'MarkerEdgeColor', C11, 'MarkerFaceColor', C11, 'MarkerSize',ms2);
plot(2, mean(AICc_h1(:,2)), 's', 'MarkerEdgeColor', C22, 'MarkerFaceColor', C22, 'MarkerSize',ms2);
plot(3, mean(AICc_h2(:,2)), 's', 'MarkerEdgeColor', C33, 'MarkerFaceColor', C33, 'MarkerSize',ms2);
plot(4, mean(AICc_h3(:,2)), 's', 'MarkerEdgeColor', C44, 'MarkerFaceColor', C44, 'MarkerSize',ms2);
errorbar(1, mean(AICc_h0(:,2)), std(AICc_h0(:,2)), 'o', 'color', C11, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
errorbar(2, mean(AICc_h1(:,2)), std(AICc_h1(:,2)), 'o', 'color', C22, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
errorbar(3, mean(AICc_h2(:,2)), std(AICc_h2(:,2)), 'o', 'color', C33, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
errorbar(4, mean(AICc_h3(:,2)), std(AICc_h3(:,2)), 'o', 'color', C44, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
xlim([0.25 4.75]); xticks(1:1:4);
xticklabels({'H0(y~x)', 'H1(y~a*x+0)', 'H2(y~1*x+b)', 'H3(y~a*x+b)'});
ylabel('AICc', 'FontName', 'Arial', 'FontSize', 10);
title('5 dots');
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 14, 'linewidth', 1.5,'TickLength',[0.025 0]);

figure
subplot(1,2,1)
plot(1, Slope(:,1), 'o', 'MarkerEdgeColor', C2, 'MarkerFaceColor', C2, 'MarkerSize',ms); hold on
plot(1, mean(Slope(:,1)), 's', 'MarkerEdgeColor', C22, 'MarkerFaceColor', C22, 'MarkerSize',ms2);
errorbar(1, mean(Slope(:,1)), std(Slope(:,1)), 'o', 'color', C22, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
plot(2, alpha(:,1), 'o', 'MarkerEdgeColor', C4, 'MarkerFaceColor', C4, 'MarkerSize',ms); hold on
plot(2, mean(alpha(:,1)), 's', 'MarkerEdgeColor', C44, 'MarkerFaceColor', C44, 'MarkerSize',ms2);
errorbar(2, mean(alpha(:,1)), std(alpha(:,1)), 'o', 'color', C44, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);

plot(3, Slope(:,2), 'o', 'MarkerEdgeColor', C2, 'MarkerFaceColor', C2, 'MarkerSize',ms); hold on
plot(3, mean(Slope(:,2)), 's', 'MarkerEdgeColor', C22, 'MarkerFaceColor', C22, 'MarkerSize',ms2);
errorbar(3, mean(Slope(:,2)), std(Slope(:,2)), 'o', 'color', C22, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
plot(4, alpha(:,2), 'o', 'MarkerEdgeColor', C4, 'MarkerFaceColor', C4, 'MarkerSize',ms); hold on
plot(4, mean(alpha(:,2)), 's', 'MarkerEdgeColor', C44, 'MarkerFaceColor', C44, 'MarkerSize',ms2);
errorbar(4, mean(alpha(:,2)), std(alpha(:,2)), 'o', 'color', C44, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);

xlim([0.25 4.75]); xticks(1:1:4);
xticklabels({'30 dots H1', '30 dots H3', '5 dots H1', '5 dots H3'});
ylabel('Estimated slope', 'FontName', 'Arial', 'FontSize', 10);
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 14, 'linewidth', 1.5,'TickLength',[0.025 0]);
lineplot(1,'h', 'k--'); 

subplot(1,2,2)
plot(1, Intercept(:,1), 'o', 'MarkerEdgeColor', C2, 'MarkerFaceColor', C2, 'MarkerSize',ms); hold on
plot(1, mean(Intercept(:,1)), 's', 'MarkerEdgeColor', C22, 'MarkerFaceColor', C22, 'MarkerSize',ms2);
errorbar(1, mean(Intercept(:,1)), std(Intercept(:,1)), 'o', 'color', C22, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
plot(2, beta(:,1), 'o', 'MarkerEdgeColor', C4, 'MarkerFaceColor', C4, 'MarkerSize',ms); hold on
plot(2, mean(beta(:,1)), 's', 'MarkerEdgeColor', C44, 'MarkerFaceColor', C44, 'MarkerSize',ms2);
errorbar(2, mean(beta(:,1)), std(beta(:,1)), 'o', 'color', C44, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);

plot(3, Intercept(:,2), 'o', 'MarkerEdgeColor', C2, 'MarkerFaceColor', C2, 'MarkerSize',ms); hold on
plot(3, mean(Intercept(:,2)), 's', 'MarkerEdgeColor', C22, 'MarkerFaceColor', C22, 'MarkerSize',ms2);
errorbar(3, mean(Intercept(:,2)), std(Intercept(:,2)), 'o', 'color', C22, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
plot(4, beta(:,2), 'o', 'MarkerEdgeColor', C4, 'MarkerFaceColor', C4, 'MarkerSize',ms); hold on
plot(4, mean(beta(:,2)), 's', 'MarkerEdgeColor', C44, 'MarkerFaceColor', C44, 'MarkerSize',ms2);
errorbar(4, mean(beta(:,2)), std(beta(:,2)), 'o', 'color', C44, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);

xlim([0.25 4.75]); xticks(1:1:4);
xticklabels({'30 dots H2', '30 dots H3', '5 dots H2', '5 dots H3'});
ylabel('Estimated intercept', 'FontName', 'Arial', 'FontSize', 10);
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 14, 'linewidth', 1.5,'TickLength',[0.025 0]);
lineplot(0,'h', 'k--'); 



