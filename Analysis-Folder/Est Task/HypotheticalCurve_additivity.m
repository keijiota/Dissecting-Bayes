clear all; close all;


x = 0:0.01:1; 

% H0: P[S]_est = P[S] + e
a = 1 ; b = 0 ;
y1 = x.*a + b ;
subplot(1,3,1)
plot(x,y1, 'k-','linewidth',1); hold on 
plot(0:0.1:1.1, 0:0.1:1.1, 'k:');
axis('square');
xlim([0 1.1]); ylim([0 1.1]); 
xticks(0:0.2:1); yticks(0:0.2:1);
xlabel('Pest[S]'); ylabel('Pest[SU]+Pest[SL]'); 
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);

% H1: P[SU] = P[SL] + b + e (b > 0)
a = 1 ;
for b = [0.1 0.25 0.4]
  y2 = x.*a + b ;
  subplot(1,3,2)
  plot(x,y2, 'k-'); hold on
end
plot(0:0.1:1.1, 0:0.1:1.1, 'k:');
axis('square');
xlim([0 1.1]); ylim([0 1.1]); 
xticks(0:0.2:1); yticks(0:0.2:1);
xlabel('Pest[S]'); ylabel('Pest[SU]+Pest[SL]'); 
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);

% H2: P[SU] = P[SL] + b + e (b < 0)
a = 1 ;
for b = [-0.1 -0.25 -0.4]
  y3 = x.*a + b ;
  
  subplot(1,3,3)
  plot(x,y3, 'k-'); hold on
end
plot(0:0.1:1.1, 0:0.1:1.1, 'k:');
axis('square');
xlim([0 1.1]); ylim([0 1.1]); 
xticks(0:0.2:1); yticks(0:0.2:1);
xlabel('Pest[S]'); ylabel('Pest[SU]+Pest[SL]'); 
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);


pos(3) = 1000; pos(4) = 400;
set(gcf, 'Position', pos);


