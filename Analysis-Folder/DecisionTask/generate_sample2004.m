clear all; close all;


x = -275:1:275; y = 1:1:3000;
[X,Y] = meshgrid(x,y);

mux = 0; muy = 0;
sdx = 75; sdy = 150; r = 0;
sdx = 50; sdy = 100; r = 0;
sdx = 100; sdy = 200; r = 0;

fxy = 1/(2*pi*sdx*sdy*sqrt(1-(r^2)))...
      * exp(-1/(1-(r^2)) * (((X-mux).^2 /2/sdx^2) + ((Y-muy).^2 /2/sdy^2) - (r*(X-mux).*(Y-muy)/sdx/sdy))) ;
  
subplot(5,4,1)
% mesh(X,Y,fxy);
% zlabel('Probability Density');
% axis('equal');
% ylim([500 2000]);  xlim([-275 275]);

view(2)
% pos(3) = 300; pos(4) = 300;
% set(gcf, 'Position', pos);

n = 30 ; 
meanmat = [mux, muy]';
covmat = [sdx^2, r; r, sdy^2];


mux_t = 20; muy_t = 20;
sdx_t = 1; sdy_t = 1;
raw_t = 0.05; 

Mu30 = []; Cov30 = [];
for i = 1:10
isstop = 0; 
while isstop == 0
s30 = ndRandn(meanmat,covmat,30); 
s5 = ndRandn(meanmat,covmat,5);

mux30 = mean(s30(1,:)); muy30 = mean(s30(2,:));
covxy30 = cov(s30') ; 
sdx30 = covxy30(1)^(1/2); sdy30 = covxy30(4)^(1/2); rawxy30 = corr(s30'); rawxy30 = rawxy30(1,2);

mux5 = mean(s5(1,:)); muy5 = mean(s5(2,:));
covxy5 = cov(s5') ; 
sdx5 = covxy5(1)^(1/2); sdy5 = covxy5(4)^(1/2); rawxy5 = corr(s5'); rawxy5 = rawxy5(1,2);

if abs(mux30 - mux5) < mux_t
 if abs(muy30 - muy5) < muy_t
  if abs(sdx30 - sdx5) < sdx_t 
   if abs(sdy30 - sdy5) < sdy_t 
    if abs(rawxy30 - r) < raw_t  
        isstop = 1; 
        S30(:,:,i) = s30;
        S5(:,:,i) = s5;        
        
        figure(i)
        subplot(1,2,1)
        plot(s30(1,:), s30(2,:),'ko','MarkerFaceColor','w', 'MarkerSize',5);         
        axis('equal'); box off;
        ylim([-750 750]);  xlim([-275 275]);
        subplot(1,2,2)
        plot(s5(1,:), s5(2,:),'ko','MarkerFaceColor','w', 'MarkerSize',5);         
        axis('equal'); box off;
        ylim([-750 750]);  xlim([-275 275]);
        
    end
   end
  end
 end
end
end     
end

fa
s30 = S30(:,:,1); 
s30 = s30 - mean(s30,2); 
s5 = S5(:,:,1);
s5 = s5 - mean(s5,2); 

save sampleFigure3.mat s30 s5







