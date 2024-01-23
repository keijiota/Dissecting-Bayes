clear all; close all;


x = -5:.01:5; y = -5:.01:5;
[X,Y] = meshgrid(x,y);

mux = 0; muy = 0;
sdx = 0.75; sdy = 1.5; r = 0;

fxy = 1/(2*pi*sdx*sdy*sqrt(1-(r^2)))...
      * exp(-1/(1-(r^2)) * (((X-mux).^2 /2/sdx^2) + ((Y-muy).^2 /2/sdy^2) - (r*(X-mux).*(Y-muy)/sdx/sdy))) ;
  
subplot(5,4,1)
mesh(X,Y,fxy);
zlabel('Probability Density');
axis('equal');
xlim([-3 3]); ylim([-4.7 4.7]);

view(2)
% pos(3) = 300; pos(4) = 300;
% set(gcf, 'Position', pos);

n = 30 ; 
meanmat = [mux, muy]';
covmat = [sdx^2, r; r, sdy^2];


mux_t = 0.1; muy_t = 0.1;
sdx_t = 0.1; sdy_t = 0.1;
raw_t = 0.01; 

Mu = []; Cov = [];
for i = 1:19
isstop = 0; 
while isstop == 0
s = ndRandn(meanmat,covmat,n);

musx = mean(s(1,:)); musy = mean(s(2,:));
covsxy = cov(s') ; 
sdsx = covsxy(1)^(1/2); sdsy = covsxy(4)^(1/2); rawxy = covsxy(2);

if abs(musx - mux) < mux_t
 if abs(musy - muy) < muy_t
  if abs(sdsx - sdx) < sdx_t 
   if abs(sdsy - sdy) < sdy_t 
    if abs(rawxy - r) < raw_t  
        isstop = 1; 
        S(:,:,i) = s;
        Mu = [Mu; musx, musy];
        Cov = [Cov; sdsx, sdsy, rawxy];
        subplot(5,4,i+1)
        plot(s(1,:), s(2,:),'ko','MarkerFaceColor','w', 'MarkerSize',5); 
        axis('equal'); box off;
        xlim([-3 3]); ylim([-4.7 4.7]);
    end
   end
  end
 end
end
end     
end

fa
s = S(:,:,1)
save s30_4 S




