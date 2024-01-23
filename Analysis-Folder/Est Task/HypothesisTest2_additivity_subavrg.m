clear all; close all;

load('dat_percep.mat');

xx = 0:0.01:1; result = [];
colormap = [0 0 207; 0 111 255; 0  255 255; 111 255 143;
    207 255 48; 255 191 0; 255 96 0; 200 0  0; 100  0  0];
for condition = 1:2
    if condition == 1
        % 30, symmetric
        x = MeanEstProb(:,:,1) ; y1 = MeanEstProb(:,:,2); y2 = MeanEstProb(:,:,3);
        conname = '30 dots';
    else
        % 5, symmetric
        x = MeanEstProb(:,:,4) ; y1 = MeanEstProb(:,:,5); y2 = MeanEstProb(:,:,6);
        conname = '5 dots';
    end
    X = mean(x); Y = mean(y1+y2);
    %     X = reshape(X, 1,153); Y = reshape(Y, 1,153);
    
    %         model test
    %                 a = 1; b = 0;
    %                 Y = X.*a + b + (rand(1,9)*0.0001-0.00005);
    %                 a = 0.8+rand*0.4; b = rand*0.2-0.1 ; b=0;
    %                 Y = X.*a + b ;
    %                 gamma = rand+0.5 ;
    %                 Y = (X.^gamma ./ ((X.^gamma + (1-X).^gamma).^(1/gamma))); % Tversky & Kahneman 1992
    %                 Y = exp(-(-log(X)).^gamma) ;
    
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

    %         % H3: P[SU] = a*P[SL] + b + e, parameter free model
    %         clear fitparams fitparams_pre mnll mnll_pre
    %         N = length(X) ; K = 2;
    %         for rep = 1:5
    %             NLL = @(params) [N/2 * (log(2*pi) + log(mean(((params(1)*X + params(2)) - Y).^2)) + 1)];
    %             [fitparams_pre mnll_pre] = fmincon(NLL,[rand(2,1)*2],[],[],[],[],[],[]);
    %             if rep == 1
    %                 fitparams = fitparams_pre; mnll = mnll_pre;
    %             elseif rep > 1 && mnll >= mnll_pre
    %                 fitparams = fitparams_pre; mnll = mnll_pre;
    %             end
    %         end
    %         Y_hut = fitparams(1)*X + fitparams(2);
    %         Alpha3(1,condition) = fitparams(1);
    %         Beta3(1,condition) = fitparams(2);
    %         Residual_h3(1,condition) = mean((Y-Y_hut).^2) ;
    %         nll_h3 = - length(X)/2 * (log(2*pi) + log(Residual_h3(1,condition)) + 1) ;
    %         AICc_h3(1,condition) = -2*(nll_h3) + 2*K + 2*(K*(K+1))/(N-K-1);
    %         y_h3 = xx*fitparams(1) + fitparams(2);
    
    subplot(1,2,condition)
    for i = 1:9
        errorbar(mean(x(:,i)), mean([y1(:,i)+y2(:,i)]), 2*std(x(:,i))/sqrt(17), 'horizontal', 'k', 'linewidth', 0.5, 'MarkerSize', 0.1, 'CapSize',5); hold on
        errorbar(mean(x(:,i)), mean([y1(:,i)+y2(:,i)]), [2*std(y1(:,i))/sqrt(17) + 2*std(y2(:,i))/sqrt(17)]/2, 'k', 'linewidth', 0.5, 'MarkerSize', 0.1, 'CapSize',5); hold on
        plot(mean(x(:,i)), mean([y1(:,i)+y2(:,i)]), 'ko', 'MarkerFaceColor',[colormap(i,:)/255], 'Markersize',7) ; hold on
    end
    %         plot(xx, y_h0, 'r-');
    %         plot(xx, y_h1, 'b-');
    %         plot(xx, y_h2, 'g-');
    %         plot(xx, y_h3, 'c-');
    %         plot(xx, y_h4, 'y-');
    %         plot(xx, y_h5, 'm-');
    plot(0:0.1:1.2,0:0.1:1.2,'k--');
    plot(xx, y_h1,'k-');
    lineplot(1,'h', 'k:'); lineplot(1,'v', 'k:');
    axis('square');
    xlim([0 1.05]); ylim([0 1.05]);
    xticks(0:0.2:1.2); yticks(0:0.2:1.2);
    xlabel('\pi(P[S])', 'FontName', 'Arial', 'FontSize', 10);
    ylabel('\pi(P[SU])+\pi(P[SL])', 'FontName', 'Arial', 'FontSize', 10);
    title(strcat('Sub', num2str(1),',',conname));
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
mean(Meandiff)

beta = [Beta; Beta2]

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
xlim([0.25 3.75]); xticks(1:1:6);
xticklabels({'H0(y~x)', 'H1(y~a*x+0)', 'H2(y~1*x+b)', 'H3(y~a*x+b)', 'H4(y~x^g/(x^g+(1-x)^g)^1^/^g', 'H5(y~exp(-(-ln(x)^g))'});
ylabel('AICc', 'FontName', 'Arial', 'FontSize', 10);
title('30 dots');
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 14, 'linewidth', 1.5,'TickLength',[0.025 0]);

subplot(1,2,2)
plot(1, AICc_h0(:,2), 's', 'MarkerEdgeColor', C11, 'MarkerFaceColor', C11, 'MarkerSize',ms2); hold on
plot(2, AICc_h1(:,2), 's', 'MarkerEdgeColor', C22, 'MarkerFaceColor', C22, 'MarkerSize',ms2);
plot(3, AICc_h2(:,2), 's', 'MarkerEdgeColor', C33, 'MarkerFaceColor', C33, 'MarkerSize',ms2);
xlim([0.25 3.75]); xticks(1:1:6);
xticklabels({'H0(y~x)', 'H1(y~a*x+0)', 'H2(y~1*x+b)', 'H3(y~a*x+b)', 'H4(y~x^g/(x^g+(1-x)^g)^1^/^g', 'H5(y~exp(-(-ln(x)^g))'});
ylabel('AICc', 'FontName', 'Arial', 'FontSize', 10);
title('5 dots');
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 14, 'linewidth', 1.5,'TickLength',[0.025 0]);

