clear all; close all;

%% plot time series of endpoint

load('dat_decision');
for condition = 1:4
    score = Score(:,:,condition);
%     optscore = OptEG(:,:,condition);    
    optscore = OptEEG(:,:,condition);    
    
    meanscore = nanmean(score') ; 
    meanoptscore = nanmean(optscore') ; 
    efficiency = meanscore./meanoptscore;        
    
    [N T] = size(score);    
    % mean across trials
    Efficiency(:,condition) = efficiency' ;
    miss = score<=0;    
    summiss(:, condition) = sum(miss');    
end

% bar plot 
tmp = [1 3 2 4]; % [30,0], [5,0], [30,-500], [5,-500]
Efficiency = Efficiency(:,tmp) ; 
summiss = summiss(:,tmp) ; 

ms = 9; cs = 10; ms2 = 10;
figure
C1 = [204 247 255]./255 ; C11 = [80 180 255]./255 ; 
C3 = [179 179 255]./255; C33 = [25 25 255]./255;
C2 = [255 221 200]./255; C22 = [255 140 0]./255;
C4 = [255 200 200]./255; C44 = [200 0 0]./255; 

ic = [1 2.2 3.4 4.6]; p = 0.3;

subplot(1,2,1)
plot(ic(1), summiss(:,1), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', C11, 'MarkerSize',ms); hold on
plot(ic(2), summiss(:,2), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', C22, 'MarkerSize',ms); 
plot(ic(3), summiss(:,3), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', C33, 'MarkerSize',ms); 
plot(ic(4), summiss(:,4), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', C44, 'MarkerSize',ms); 
plot(ic(1)+p, mean(summiss(:,1)), 's', 'MarkerEdgeColor', C11, 'MarkerFaceColor', C11, 'MarkerSize',ms2); hold on
plot(ic(2)+p, mean(summiss(:,2)), 's', 'MarkerEdgeColor', C22, 'MarkerFaceColor', C22, 'MarkerSize',ms2); 
plot(ic(3)+p, mean(summiss(:,3)), 's', 'MarkerEdgeColor', C33, 'MarkerFaceColor', C33, 'MarkerSize',ms2); 
plot(ic(4)+p, mean(summiss(:,4)), 's', 'MarkerEdgeColor', C44, 'MarkerFaceColor', C44, 'MarkerSize',ms2); 
errorbar(ic(1)+p, mean(summiss(:,1)), 2*std(summiss(:,1))/sqrt(17), 'o', 'color', C11, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs); 
errorbar(ic(2)+p, mean(summiss(:,2)), 2*std(summiss(:,2))/sqrt(17), 'o', 'color', C22, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs); 
errorbar(ic(3)+p, mean(summiss(:,3)), 2*std(summiss(:,3))/sqrt(17), 'o', 'color', C33, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs); 
errorbar(ic(4)+p, mean(summiss(:,4)), 2*std(summiss(:,4))/sqrt(17), 'o', 'color', C44, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs); 
plot(ic(1:2)+p, [mean(summiss(:,1)), mean(summiss(:,2))], 'k--', 'linewidth', 1); 
plot(ic(3:4)+p, [mean(summiss(:,3)), mean(summiss(:,4))], 'k--', 'linewidth', 1);   
ylim([0 20]);  xlim([ic(1)-0.75 ic(4)+0.75]); xticks(ic+p/2); yticks(0:4:20); 
xticklabels({'30 dots, 0 penalty', '5 dots, 0 penalty', '30 dots, -500 penalty', '5 dots, -500 penalty'}); 
ylabel('Numbers of mistrial', 'FontName', 'Arial', 'FontSize', 10);
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);

subplot(1,2,2)
plot(ic(1), Efficiency(:,1), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', C11, 'MarkerSize',ms); hold on
plot(ic(2), Efficiency(:,2), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', C22, 'MarkerSize',ms); 
plot(ic(3), Efficiency(:,3), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', C33, 'MarkerSize',ms); 
plot(ic(4), Efficiency(:,4), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', C44, 'MarkerSize',ms); 
plot(ic(1)+p, mean(Efficiency(:,1)), 's', 'MarkerEdgeColor', C11, 'MarkerFaceColor', C11, 'MarkerSize',ms2); hold on
plot(ic(2)+p, mean(Efficiency(:,2)), 's', 'MarkerEdgeColor', C22, 'MarkerFaceColor', C22, 'MarkerSize',ms2); 
plot(ic(3)+p, mean(Efficiency(:,3)), 's', 'MarkerEdgeColor', C33, 'MarkerFaceColor', C33, 'MarkerSize',ms2); 
plot(ic(4)+p, mean(Efficiency(:,4)), 's', 'MarkerEdgeColor', C44, 'MarkerFaceColor', C44, 'MarkerSize',ms2); 
errorbar(ic(1)+p, mean(Efficiency(:,1)), 2*std(Efficiency(:,1))/sqrt(17), 'o', 'color', C11, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs); 
errorbar(ic(2)+p, mean(Efficiency(:,2)), 2*std(Efficiency(:,2))/sqrt(17), 'o', 'color', C22, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs); 
errorbar(ic(3)+p, mean(Efficiency(:,3)), 2*std(Efficiency(:,3))/sqrt(17), 'o', 'color', C33, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs); 
errorbar(ic(4)+p, mean(Efficiency(:,4)), 2*std(Efficiency(:,4))/sqrt(17), 'o', 'color', C44, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs); 
plot(ic(1:2)+p, [mean(Efficiency(:,1)), mean(Efficiency(:,2))], 'k--', 'linewidth', 1); 
plot(ic(3:4)+p, [mean(Efficiency(:,3)), mean(Efficiency(:,4))], 'k--', 'linewidth', 1);   
lineplot(1, 'h', 'k-', 'linewidth',1); 
ylim([-0.45 1.2]);xlim([ic(1)-0.75 ic(4)+0.75]); xticks(ic+p/2); yticks(-0.4:0.2:1.2); 
xticklabels({'30 dots, 0 penalty', '5 dots, 0 penalty', '30 dots, -500 penalty', '5 dots, -500 penalty'}); 
ylabel('Efficiency', 'FontName', 'Arial', 'FontSize', 10);
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);

pos(3) = 1050; pos(4) = 450;
set(gcf, 'Position', pos);

Result = [];
for i = 1:4
    [tvalue, df, p, h, d, meandiff, CI_l, CI_u] = myttest(Efficiency(:,i)*100,100);
    h = p < 0.05/4;
    Result = [Result; tvalue, df, p, h, d, meandiff, CI_l, CI_u];    
end
[tvalue, df, p, h, d, meandiff, CI_l, CI_u] = myttest(Efficiency(:,1)*100, Efficiency(:,2)*100);
Result = [Result; tvalue, df, p, h, d, meandiff, CI_l, CI_u];
[tvalue, df, p, h, d, meandiff, CI_l, CI_u] = myttest(Efficiency(:,3)*100, Efficiency(:,4)*100);
Result = [Result; tvalue, df, p, h, d, meandiff, CI_l, CI_u];
mean(Efficiency(:,:)*100)
mean(mean(Efficiency(:,:)*100))

Result2 = [];
[tvalue, df, p, h, d, meandiff, CI_l, CI_u] = myttest(summiss(:,2), summiss(:,1));
Result2 = [Result2; tvalue, df, p, h, d, meandiff, CI_l, CI_u];
[tvalue, df, p, h, d, meandiff, CI_l, CI_u] = myttest(summiss(:,4), summiss(:,3));
Result2 = [Result2; tvalue, df, p, h, d, meandiff, CI_l, CI_u];



