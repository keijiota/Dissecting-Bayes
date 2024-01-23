clear all; close all;

%% plot time series of endpoint

load('dat_decision');
for condition = 1:4
    endpointorg = EndpointorgY(:,:,condition); 
    
    [N T] = size(endpointorg);        
    
    meanendpoint(:,condition) = nanmean(endpointorg) ;
end


tmp = [1 3 2 4]; % [30,0], [5,0], [30,-500], [5,-500]
meanendpoint = meanendpoint(:,tmp) ; 


ms = 6.5; cs = 10; ms2 = 10;
C1 = [204 247 255]./255 ; C11 = [80 180 255]./255 ; 
C3 = [179 179 255]./255; C33 = [255 140 0]./255;
C2 = [255 221 200]./255; C22 = [25 25 255]./255;
C4 = [255 200 200]./255; C44 = [200 0 0]./255; 

figure
subplot(2,2,1)
plot(1:50, meanendpoint(:,1), 'o', 'MarkerEdgeColor', C11, 'MarkerFaceColor', C11, 'MarkerSize',ms); hold on
lineplot(mean(meanendpoint(:,1)), 'h','--', 'Color', 'k', 'linewidth',1); 
ylim([120 170]); yticks(100:10:180); xlim([0 51]); xticks(0:10:50);
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10, 'linewidth', 1,'TickLength',[0.012 0]);

subplot(2,2,2)
plot(1:50, meanendpoint(:,2), 'o', 'MarkerEdgeColor', C22, 'MarkerFaceColor', C22, 'MarkerSize',ms); hold on
lineplot(mean(meanendpoint(:,2)), 'h','--', 'Color', 'k',  'linewidth',1); 
ylim([120 170]); yticks(100:10:180); xlim([0 51]); xticks(0:10:50);
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10, 'linewidth', 1,'TickLength',[0.012 0]);

subplot(2,2,3)
plot(1:50, meanendpoint(:,3), 'o', 'MarkerEdgeColor', C33, 'MarkerFaceColor', C33, 'MarkerSize',ms); hold on
lineplot(mean(meanendpoint(:,3)), 'h','--', 'Color', 'k', 'linewidth',1); 
ylim([120 170]); yticks(100:10:180); xlim([0 51]); xticks(0:10:50);
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10, 'linewidth', 1,'TickLength',[0.012 0]);

subplot(2,2,4)
plot(1:50, meanendpoint(:,4), 'o', 'MarkerEdgeColor', C44, 'MarkerFaceColor', C44, 'MarkerSize',ms); hold on
lineplot(mean(meanendpoint(:,4)), 'h','--', 'Color', 'k',  'linewidth',1); 
ylim([120 170]); yticks(100:10:180); xlim([0 51]); xticks(0:10:50);
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10, 'linewidth', 1,'TickLength',[0.012 0]);


pos(3) = 1550; pos(4) = 500;
set(gcf, 'Position', pos);


