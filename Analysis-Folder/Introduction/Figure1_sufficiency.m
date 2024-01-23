clear all; close all;


n = 30 ; 
meanmat = [mux, muy]';
covmat = [sdx^2, r; r, sdy^2];

% mux_t = 0.2; muy_t = 0.2;
% sdx_t = 1.5; sdy_t = 3;
% raw_t = 0.15; 
% isstop = 0; 
% while isstop == 0
% s = ndRandn(meanmat,covmat,n);
% 
% musx = mean(s(1,:)); musy = mean(s(2,:));
% covsxy = cov(s') ; 
% sdsx = covsxy(1)^(1/2); sdsy = covsxy(4)^(1/2); rawxy = covsxy(2);
% 
% if abs(musx - mux) < mux_t
%  if abs(musy - muy) < muy_t
%   if abs(sdsx - sdx) < sdx_t 
%    if abs(sdsy - sdy) < sdy_t 
%     if abs(rawxy - r) < raw_t  
%         isstop = 1; 
%     end
%    end
%   end
%  end
% end
% end     

load('s30_4'); 
s = S(:,:,10);
musx = mean(s(1,:)); musy = mean(s(2,:));
covsxy = cov(s') ; 
sdsx = covsxy(1)^(1/2); sdsy = covsxy(4)^(1/2); rawxy = covsxy(2);

figure
plot(s(1,:), s(2,:),'ko','MarkerFaceColor','w', 'MarkerSize',10); 
axis('equal');
xlim([-2.7 2.7]); ylim([-4.5 4.5]); 
xticks(-10:20:10); yticks(-10:20:10);
pos(3) = 550; pos(4) = 550;
set(gcf, 'Position', pos);

fxyh = 1/(2*pi*sdsx*sdsy*sqrt(1-(rawxy^2)))...
      * exp(-1/(1-(rawxy^2)) * (((X-musx).^2 /2/sdsx^2) + ((Y-musy).^2 /2/sdsy^2) - (rawxy*(X-musx).*(Y-musy)/sdsx/sdsy)));  

figure
mesh(X,Y,fxyh);
zlabel('Probability Density');
axis('equal');
xlim([-2.7 2.7]); ylim([-4.5 4.5]); 
xticks(-10:20:10); yticks(-10:20:10);
view(2);  
pos(3) = 550; pos(4) = 550;
set(gcf, 'Position', pos);


fa
save s30_3 s


  