clear all; close all;

x = -5:.01:5; y = -5:.01:5;
[X,Y] = meshgrid(x,y);

mux = 0; muy = 0;
sdx = 0.75; sdy = 1.5; r = 0;

fxy = 1/(2*pi*sdx*sdy*sqrt(1-(r^2)))...
      * exp(-1/(1-(r^2)) * (((X-mux).^2 /2/sdx^2) + ((Y-muy).^2 /2/sdy^2) - (r*(X-mux).*(Y-muy)/sdx/sdy))) ;
  
% subplot(1,4,1)
% mesh(X,Y,fxy);
% zlabel('Probability Density');
% axis('equal');
% % axis([-3 3 -4.7 4.7 0 .5]); 
% xlim([-3 3]); ylim([-4.7 4.7]);
% view(2);

s1 = load('s1.mat'); s1 = s1.s; 
s2 = load('s2_2.mat'); s2 = s2.s; 
s4 = load('s3_3.mat'); s4 = s4.s; 

s2(1,2) = s2(1,2)-0.15;
s2(1,8) = s2(1,8)+0.15;
s3 = [s2(1,:); -s2(2,:)];

s4(2,5) = s4(2,5)+0.2;
s4(2,2) = s4(2,2)-0.2;


[mean(s1,2) mean(s2,2) mean(s3,2) mean(s4,2)]
[cov(s1'),cov(s2'),cov(s3'),cov(s4')]


subplot(1,5,2)
plot(s1(1,:), s1(2,:),'ko','MarkerFaceColor','w', 'MarkerSize',8);
axis('equal'); box off;
xlim([-3 3]); ylim([-4.7 4.7]);

subplot(1,5,3)
plot(s2(1,:), s2(2,:),'ko','MarkerFaceColor','w', 'MarkerSize',8);
axis('equal'); box off;
xlim([-3 3]); ylim([-4.7 4.7]);

subplot(1,5,4)
plot(s3(1,:), s3(2,:),'ko','MarkerFaceColor','w', 'MarkerSize',8);
axis('equal'); box off;
xlim([-3 3]); ylim([-4.7 4.7]);

subplot(1,5,5)
plot(s4(1,:), s4(2,:),'ko','MarkerFaceColor','w', 'MarkerSize',8);
axis('equal'); box off;
xlim([-3 3]); ylim([-4.7 4.7]);


pos(3) = 1300; pos(4) = 1300;
set(gcf, 'Position', pos);




