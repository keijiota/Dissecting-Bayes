clear; close all;

% setting
y = -500:0.1:750;
v_hut = [0:1:7500];
Y = 50:0.1:200;

TL = 180 ;
for i = 1:length(y)
    if y(i) < 0
        gainFunc(i) = 0;
        gainFunc2(i) = 0;
    elseif y(i) >= 0 && y(i) <= TL
        gainFunc(i) = y(i)/(TL)*100;
        gainFunc2(i) = y(i)/(TL)*100;
    elseif y(i) > TL
        gainFunc(i) = 0;
        gainFunc2(i) = -500;
    end
end
gainfunc = gainFunc;

load sampleFigure1
dat30 = s30'; dat5 = s5';
dat30 = dat30*0.91; dat5 = dat5*0.91;
% dat30 = dat30*0.95; dat5 = dat5*0.95;


[mean(dat30), mean(dat5)]
[cov(dat30), cov(dat5)]


ms = 5; setpoint = 135; datfb = [-3 160]; ls = 0.5;
figure
subplot(1,3,1)
plot(dat30(:,1), dat30(:,2), 'ok', 'MarkerFaceColor','w','MarkerSize',ms,'linewidth',ls); hold on
plot(0,0, 'ko', 'MarkerFaceColor','b','MarkerSize',ms+1,'linewidth',1); hold on
xlim([-80, 30]); ylim([-40 220]); xticks(-200:300:100); yticks(-200:500:300);
lineplot(180, 'h', '-g','linewidth',2);
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10, 'linewidth', 1);

subplot(1,3,2)
plot(0,0, 'ko', 'MarkerFaceColor','b','MarkerSize',ms+1,'linewidth',1); hold on
plot(dat30(:,1), dat30(:,2)+setpoint, 'ok', 'MarkerFaceColor','w','MarkerSize',ms,'linewidth',ls);
xlim([-80, 30]); ylim([-40 220]); xticks(-200:300:100); yticks(-200:500:300);
lineplot(180, 'h', '-g','linewidth',2);
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10, 'linewidth', 1);

subplot(1,3,3)
plot(0,0, 'ko', 'MarkerFaceColor','b','MarkerSize',ms+1,'linewidth',1); hold on
plot(dat30(:,1), dat30(:,2)+setpoint, 'ok', 'MarkerFaceColor','w','MarkerSize',ms,'linewidth',ls);
plot(datfb(1), datfb(2), 'ok', 'MarkerFaceColor','y','MarkerSize',ms+1.5,'linewidth',ls);
xlim([-80, 30]); ylim([-40 220]); xticks(-200:300:100); yticks(-200:500:300);
lineplot(180, 'h', '-g','linewidth',2);
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 10, 'linewidth', 1);

pos(3) = 800; pos(4) = 370;
set(gcf, 'Position', pos);

std(dat30)

xx = -50:0.1:50; yy = y; [XX,YY] = meshgrid(xx,yy);
c = jet; c = c(1:200,:);

mu1 = 0; mu2 = setpoint;
s1 = 10; s2 = 20; r = 0 ;
fxy = 1/(2*pi*s1*s2*sqrt(1-(r^2)))...
    * exp(-1/(1-(r^2)) * (((XX-mu1).^2 /2/s1^2) + ((YY-mu2).^2 /2/s2^2) - (r*(XX-mu1).*(YY-mu2)/s1/s2))) ;

figure
lineplot(180, 'h', '-g','linewidth',1);
mesh(XX,YY,fxy); hold on
axis('equal'); colormap(c); 

view(2);
xlim([-30, 30]); ylim([-40 220]); xticks(-200:300:100); yticks(-200:500:300);

pos(3) = 800; pos(4) = 370;
set(gcf, 'Position', pos);

print(gcf,'pdf.png','-dpng','-r600');       


fa
figure; ms = 4.5;
scale1 = 7; scale2 = 10;
load sampleFigure2
subplot(2,2,1);
s30 = s30 ./scale1;
plot(s30(1,:),s30(2,:), 'ok', 'MarkerFaceColor','w','MarkerSize',ms,'linewidth',ls);
xlim([-80, 30]); ylim([-70 220]); xticks(-200:300:100); yticks(-200:500:300);

subplot(2,2,2);
s5 = s5 ./scale1;
plot(s5(1,:),s5(2,:), 'ok', 'MarkerFaceColor','w','MarkerSize',ms,'linewidth',ls);
xlim([-80, 30]); ylim([-70 220]); xticks(-200:300:100); yticks(-200:500:300);

load sampleFigure3
subplot(2,2,3);
s30 = s30 ./scale2;
plot(s30(1,:),s30(2,:), 'ok', 'MarkerFaceColor','w','MarkerSize',ms,'linewidth',ls);
xlim([-80, 30]); ylim([-70 220]); xticks(-200:300:100); yticks(-200:500:300);

subplot(2,2,4);
s5 = s5 ./scale2;
plot(s5(1,:),s5(2,:), 'ok', 'MarkerFaceColor','w','MarkerSize',ms,'linewidth',ls);
xlim([-80, 30]); ylim([-70 220]); xticks(-200:300:100); yticks(-200:500:300);

pos(3) = 700; pos(4) = 770;
set(gcf, 'Position', pos);

