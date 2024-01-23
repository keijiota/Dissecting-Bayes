clear all; close all;

load('dat_percep.mat');

xx = 0:0.01:1;
colormap = [0 0 207; 0 111 255; 0  255 255; 111 255 143;
            207 255 48; 255 191 0; 255 96 0; 200 0  0; 100  0  0]; 
for condition = 1:2
    for subi = 1:17
        if condition == 1
            % 30, symmetric
            X = MeanEstProb(subi,:,1) ; Y = MeanEstProb(subi,:,2) + MeanEstProb(subi,:,3);
            figure(1)
            conname = '30 dots';
        else
            % 5, symmetric
            X = MeanEstProb(subi,:,4) ; Y = MeanEstProb(subi,:,5) + MeanEstProb(subi,:,6);
            figure(2)
            conname = '5 dots';
        end
        
        %         model test
        %                 a = 1; b = 0;
        %                 Y = X.*a + b + (rand(1,9)*0.0001-0.00005);
        %                 a = 0.8+rand*0.4; b = rand*0.2-0.1 ; b=0;
        %                 Y = X.*a + b ;
        
        % H0: P[S]_est = P[S] + e
        N = length(X) ; K = 0;
        a = 1 ; b = 0 ;
        Y_hut = X.*a + b ;
        Residual_h0(subi,condition) = mean((Y-Y_hut).^2) ;
        nll_h0 = - length(X)/2 * (log(2*pi) + log(Residual_h0(subi,condition)) + 1) ;% bigger is good
        AICc_h0(subi,condition) = -2*(nll_h0) + 2*K + 2*(K*(K+1))/(N-K-1) ; % smaller is good
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
        Alpha(subi,condition) = fitparams(1);
        Beta(subi,condition) = fitparams(2);        
        Residual_h1(subi,condition) = mean((Y-Y_hut).^2) ;
        nll_h1 = - length(X)/2 * (log(2*pi) + log(Residual_h1(subi,condition)) + 1) ;
        AICc_h1(subi,condition) = -2*(nll_h1) + 2*K + 2*(K*(K+1))/(N-K-1);
        y_h1 = xx*fitparams(1) + fitparams(2);
        
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
        Alpha2(subi,condition) = fitparams(1);
        Beta2(subi,condition) = fitparams(2);        
        Residual_h2(subi,condition) = mean((Y-Y_hut).^2) ;
        nll_h2 = - length(X)/2 * (log(2*pi) + log(Residual_h2(subi,condition)) + 1) ;
        AICc_h2(subi,condition) = -2*(nll_h2) + 2*K + 2*(K*(K+1))/(N-K-1);
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
%         Alpha3(subi,condition) = fitparams(1);
%         Beta3(subi,condition) = fitparams(2);        
%         Residual_h3(subi,condition) = mean((Y-Y_hut).^2) ;
%         nll_h3 = - length(X)/2 * (log(2*pi) + log(Residual_h3(subi,condition)) + 1) ;
%         AICc_h3(subi,condition) = -2*(nll_h3) + 2*K + 2*(K*(K+1))/(N-K-1);
%         y_h3 = xx*fitparams(1) + fitparams(2);        

        subplot(3,6,subi)
        for i = 1:9
            plot(X(i),Y(i), 'ko', 'MarkerFaceColor',[colormap(i,:)/255], 'Markersize',10) ; hold on
        end        
        %         plot(xx, y_h0, 'r-');
        %         plot(xx, y_h1, 'b-');
        %         plot(xx, y_h2, 'g-');
        %         plot(xx, y_h3, 'c-');
        %         plot(xx, y_h4, 'y-');
        %         plot(xx, y_h5, 'm-');
        plot(0:0.1:1.3,0:0.1:1.3,'k--');
        lineplot(1,'h','k:'); lineplot(1,'v','k:');        
        axis('square');
        xlim([0 1.3]); ylim([0 1.3]);
        xticks(0:0.5:1.3); yticks(0:0.5:1.3);
%         xlabel('\pi(P[S])', 'FontName', 'Arial', 'FontSize', 10);
%         ylabel('\pi(P[SU])+\pi(P[SL])', 'FontName', 'Arial', 'FontSize', 10);
        title(strcat('Sub', num2str(subi),',',conname));
        set(gca, 'Fontname', 'Arial', 'Fontsize', 14, 'linewidth', 1, 'TickLength',[0.025 0]);
    end
end

deltaAicc2null_30d = AICc_h0(:,1) - [AICc_h0(:,1) AICc_h1(:,1) AICc_h2(:,1)];
% negative means null is better
evidenceratio_30d = exp((deltaAicc2null_30d - deltaAicc2null_30d(1))/2) ;
% er of 0.33 means null is 3 times likely, er of 3 means null is 3 times less likely
deltaAicc2null_5d = AICc_h0(:,2) - [AICc_h0(:,2) AICc_h1(:,2) AICc_h2(:,2)];
evidenceratio_5d = exp((deltaAicc2null_5d - deltaAicc2null_5d(1))/2) ;

meandeltaAicc2null_30d =  mean(AICc_h0(:,1)) - [mean(AICc_h0(:,1)) mean(AICc_h1(:,1)) mean(AICc_h2(:,1))];
meanevidenceratio_30d = exp((meandeltaAicc2null_30d - meandeltaAicc2null_30d(1))/2) ;
meandeltaAicc2null_5d = mean(AICc_h0(:,2)) - [mean(AICc_h0(:,2)) mean(AICc_h1(:,2)) mean(AICc_h2(:,2))];
meanevidenceratio_5d = exp((meandeltaAicc2null_5d - meandeltaAicc2null_5d(1))/2) ;

result = [meandeltaAicc2null_30d', meandeltaAicc2null_5d', meanevidenceratio_30d', meanevidenceratio_5d']; 

result2 = [mean(Alpha) mean(Alpha2); mean(Beta) mean(Beta2)];


% fa
% C1 = [255 200 200]./255; C11 = [200 0 0]./255;
% C2 = [179 179 255]./255; C22 = [25 25 255]./255;
% C3 = [200 255 200]./255 ; C33 = [52 255 50]./255 ;
% C4 = [204 247 255]./255 ; C44 = [80 180 255]./255 ;
% C5 = [255 240 0]./255 ; C55 = [255 155 0]./255 ;
% C6 = [240 210 240]./255 ; C66 = [255 0 255]./255 ;
% ms = 8; cs = 10; ms2 = 10;
% 
% figure
% subplot(1,2,1)
% plot(1, AICc_h0(:,1), 'o', 'MarkerEdgeColor', C1, 'MarkerFaceColor', C1, 'MarkerSize',ms); hold on
% plot(2, AICc_h1(:,1), 'o', 'MarkerEdgeColor', C2, 'MarkerFaceColor', C2, 'MarkerSize',ms);
% plot(3, AICc_h2(:,1), 'o', 'MarkerEdgeColor', C3, 'MarkerFaceColor', C3, 'MarkerSize',ms);
% plot(1, mean(AICc_h0(:,1)), 's', 'MarkerEdgeColor', C11, 'MarkerFaceColor', C11, 'MarkerSize',ms2);
% plot(2, mean(AICc_h1(:,1)), 's', 'MarkerEdgeColor', C22, 'MarkerFaceColor', C22, 'MarkerSize',ms2);
% plot(3, mean(AICc_h2(:,1)), 's', 'MarkerEdgeColor', C33, 'MarkerFaceColor', C33, 'MarkerSize',ms2);
% errorbar(1, mean(AICc_h0(:,1)), std(AICc_h0(:,1)), 'o', 'color', C11, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
% errorbar(2, mean(AICc_h1(:,1)), std(AICc_h1(:,1)), 'o', 'color', C22, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
% errorbar(3, mean(AICc_h2(:,1)), std(AICc_h2(:,1)), 'o', 'color', C33, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
% xlim([0.25 3.75]); xticks(1:1:6);
% xticklabels({'H0(y~x)', 'H1(y~a*x+0)', 'H2(y~1*x+b)', 'H3(y~a*x+b)', 'H4(y~x^g/(x^g+(1-x)^g)^1^/^g', 'H5(y~exp(-(-ln(x)^g))'});
% ylabel('AICc', 'FontName', 'Arial', 'FontSize', 10);
% title('30 dots');
% set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 14, 'linewidth', 1.5,'TickLength',[0.025 0]);
% 
% subplot(1,2,2)
% plot(1, AICc_h0(:,2), 'o', 'MarkerEdgeColor', C1, 'MarkerFaceColor', C1, 'MarkerSize',ms); hold on
% plot(2, AICc_h1(:,2), 'o', 'MarkerEdgeColor', C2, 'MarkerFaceColor', C2, 'MarkerSize',ms);
% plot(3, AICc_h2(:,2), 'o', 'MarkerEdgeColor', C3, 'MarkerFaceColor', C3, 'MarkerSize',ms);
% plot(1, mean(AICc_h0(:,2)), 's', 'MarkerEdgeColor', C11, 'MarkerFaceColor', C11, 'MarkerSize',ms2);
% plot(2, mean(AICc_h1(:,2)), 's', 'MarkerEdgeColor', C22, 'MarkerFaceColor', C22, 'MarkerSize',ms2);
% plot(3, mean(AICc_h2(:,2)), 's', 'MarkerEdgeColor', C33, 'MarkerFaceColor', C33, 'MarkerSize',ms2);
% errorbar(1, mean(AICc_h0(:,2)), std(AICc_h0(:,2)), 'o', 'color', C11, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
% errorbar(2, mean(AICc_h1(:,2)), std(AICc_h1(:,2)), 'o', 'color', C22, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
% errorbar(3, mean(AICc_h2(:,2)), std(AICc_h2(:,2)), 'o', 'color', C33, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
% xlim([0.25 3.75]); xticks(1:1:6);
% xticklabels({'H0(y~x)', 'H1(y~a*x+0)', 'H2(y~1*x+b)', 'H3(y~a*x+b)', 'H4(y~x^g/(x^g+(1-x)^g)^1^/^g', 'H5(y~exp(-(-ln(x)^g))'});
% ylabel('AICc', 'FontName', 'Arial', 'FontSize', 10);
% title('5 dots');
% set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 14, 'linewidth', 1.5,'TickLength',[0.025 0]);
% 
