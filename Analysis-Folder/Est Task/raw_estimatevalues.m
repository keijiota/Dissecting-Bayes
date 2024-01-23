clear; close all;


load('dat_percep.mat');
colormap = [0 0 207; 0 111 255; 0  255 255; 111 255 143;
            207 255 48; 255 191 0; 255 96 0; 200 0  0; 100  0  0];

% MeanEstProb(subi,:,1) ; % S, 30       
% MeanEstProb(subi,:,2) ; % SU, 30
% MeanEstProb(subi,:,3) ; % SL, 30        
% MeanEstProb(subi,:,4) ; % S, 5
% MeanEstProb(subi,:,5) ; % SU, 5
% MeanEstProb(subi,:,6) ; % SL, 5      

%% 
close all 
subpi = 0;
conname = {'30 dots, symmetric', '5 dots, symmetric'};

figure; k = 0;

for condition = [1 4] % 30, symmetric; 5, symmetric
    estprob = sEstProb(:,:,condition);   
    % 17 subjects * 45 trials * 6 conditions 
    % columns 1-5: 0.1 probability trials 1-5
    % columns 6-10: 0.2 probability trials 1-5
    % columns 40-45: 0.9 probability trials 1-5
    for prob = [2 8] % must have only 0.2, 0.8 in 5dots
    if prob == 2
        xl = [0 0.65];
    else 
        xl = [0.35 1]; 
    end        
    k = k +1;    
    subplot(2,2,k); hold on     
    for subi = 1:17        
        trialwise_estprob = estprob(subi, prob*5-4:prob*5);
        for itr = 1:5
            plot(trialwise_estprob(1,itr), 1*subi, 'k|','MarkerSize',5, 'LineWidth', 1);            
        end
    end
    cl = [0.6 0.6 0.6];
    if condition == 4 
        xline(0, '-',  'Color',cl,'LineWidth',1);        
        xline(0.2, '-', 'Color',cl,'LineWidth',1);
        xline(0.4, '-',  'Color',cl,'LineWidth',1);
        xline(0.6, '-',  'Color',cl,'LineWidth',1);
        xline(0.8, '-',  'Color',cl,'LineWidth',1);    
        xline(1, '-',  'Color',cl,'LineWidth',1);        
    end
    xline(prob*0.1, 'r--', 'LineWidth',2);    
    yticks(0:5:17); ylim([0 18]); xlim(xl); xticks(0:0.1:1);
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10.5, 'linewidth', 1,'TickLength',[0.015 0]);
    box off;
    
    end
end


pos(3) = 1000; pos(4) = 400;
set(gcf, 'Position', pos);

