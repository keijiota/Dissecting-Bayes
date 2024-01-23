clear all; close all;

% plot repetition vs probability estimation

load('dat_percep.mat');
colormap = [0 0 207; 0 111 255; 0  255 255; 111 255 143;
            207 255 48; 255 191 0; 255 96 0; 200 0  0; 100  0  0];

subpi = 0;
conname = {'30 dots, symmetric', '5 dots, symmetric'};
for condition = [1 4] % 30, symmetric; 5, symmetric
    subpi = subpi + 1;
    estprob = sEstProb(:,:,condition);   
    
    [N T] = size(estprob);  Nrep = 5;
    % plot average
    ms = 8; 
    subplot(1,2, subpi)        
    for ii = 1:9
        plot(1:5, mean(estprob(:,ii*5-4:ii*5)), 'ko-', 'MarkerFaceColor',[colormap(ii,:)]./255,'MarkerSize',ms); hold on
%         errorbar(1:5, mean(estprob(:,ii*5-4:ii*5)), std(estprob(:,ii*5-4:ii*5))/sqrt(N), 'ko', 'MarkerFaceColor', 'w', 'linewidth', 1, 'MarkerSize', 0.1, 'CapSize',0) ; hold on
    end

    xlim([0.5 5.5]); ylim([0 1]); xticks(1:1:5);  yticks(0:0.2:1);
    xlabel('Number of repetitions', 'FontName', 'Arial', 'FontSize', 10);
    ylabel('Estimated probability', 'FontName', 'Arial', 'FontSize', 10);
    title(conname(subpi));    
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);
end

pos(3) = 1000; pos(4) = 400;
set(gcf, 'Position', pos);

