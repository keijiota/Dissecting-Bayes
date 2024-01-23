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

% H1: P[S]_est = alpha * P[S] + beta + e
b = 0 ;
for a = [0.2 0.4 0.6 1/0.6 1/0.4 1/0.2]
  y2 = x.*a + b ;
  subplot(1,4,2)
  plot(x,y2, 'k-'); hold on
end
plot(0:0.1:1.1, 0:0.1:1.1, 'k:');
axis('square');
xlim([0 1.1]); ylim([0 1.1]); 
xticks(0:0.2:1); yticks(0:0.2:1);
xlabel('P[S]'); ylabel('Pest[S]'); 
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);

% H2: P[S]_est = P[S]^gamma / (P[S]^gamma + 1-P[S]^gamma)^(1/gamma), Tversky & Kahneman 1992
for gamma = [0.4 0.6 0.8 1/0.8 1/0.6 1/0.4]
  y3 = x.^gamma./((x.^gamma + (1-x).^gamma).^(1/gamma)) ;

  subplot(1,4,3)
  plot(x,y3, 'k-'); hold on
end
plot(0:0.1:1.1, 0:0.1:1.1, 'k:');
axis('square');
xlim([0 1.1]); ylim([0 1.1]); 
xticks(0:0.2:1); yticks(0:0.2:1);
xlabel('P[S]'); ylabel('Pest[S]'); 
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);


% H3: P[S]_est = exp(-(-ln(p)^gamma)),Prelec 1998

for gamma = [0.3 0.5 0.8 1/0.8 1/0.5 1/0.3]
  y4 = exp(-(-log(x)).^gamma) ;
  
  subplot(1,4,4)
  plot(x,y4, 'k-'); hold on
end
plot(0:0.1:1.1, 0:0.1:1.1, 'k:');
axis('square');
xlim([0 1.1]); ylim([0 1.1]); 
xticks(0:0.2:1); yticks(0:0.2:1);
xlabel('P[S]'); ylabel('Pest[S]'); 
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);


pos(3) = 1400; pos(4) = 400;
set(gcf, 'Position', pos);


