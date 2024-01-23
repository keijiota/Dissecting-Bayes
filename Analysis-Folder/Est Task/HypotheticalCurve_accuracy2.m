clear all; close all;


x = 0:0.01:1; 

% H0: P[S]_est = P[S] + e
a = 1 ; b = 0 ;
y1 = x.*a + b ;
subplot(1,4,1)
plot(x,y1, 'k-','linewidth',1); hold on 
plot(0:0.1:1.1, 0:0.1:1.1, 'k:');
axis('square');
xlim([0 1.1]); ylim([0 1.1]); 
xticks(0:0.2:1); yticks(0:0.2:1);
xlabel('P[S]'); ylabel('Pest[S]'); 
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);

% H1: LLO
for gmma = [0.3 0.6 1.4 2.6]
    lop = log(x./(1-x)) ;
    p0 = 0.6;
    lop0 = log(p0./(1-p0)) ;
    lopd = gmma*lop + (1-gmma).*lop0;    
    y2 = exp(lopd) ./ (1+exp(lopd));
    subplot(1,4,2)
    plot(x,y2, 'k-'); hold on
end
plot(0:0.1:1.1, 0:0.1:1.1, 'k:');
axis('square');
xlim([0 1.1]); ylim([0 1.1]); 
xticks(0:0.2:1); yticks(0:0.2:1);
xlabel('P[S]'); ylabel('Pest[S]'); 
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);

pos(3) = 1400; pos(4) = 400;
set(gcf, 'Position', pos);