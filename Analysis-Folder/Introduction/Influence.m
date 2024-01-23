clear; close all;


n = 8 ; 
meanmat = [mux, muy]';
covmat = [sdx^2, r; r, sdy^2];

mux_t = 0.1; muy_t = 0.1;
sdx_t = 0.1; sdy_t = 0.1;
raw_t = 0.05; 
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
    end
   end
  end
 end
end
end     

% save s8_6 s

load('s8_5'); 
mus = mean(s,2)
covsxy = cov(s');
[covsxy(1,1)^(1/2) covsxy(1,2) covsxy(2,2)^(1/2)]

sx = 0.45; sy = 0.3;

% figure
plot(s(1,:), s(2,:),'ko','MarkerFaceColor','w', 'MarkerSize',9); hold on
plot(mux, muy, 'bs','MarkerFaceColor','b', 'MarkerSize',9);
plot(sx, sy, 'rd','MarkerFaceColor','r', 'MarkerSize',9);
axis('equal');
xlim([-2.7 2.7]); ylim([-4.5 4.5]); 
xticks(-10:20:10); yticks(-10:20:10);
pos(3) = 550; pos(4) = 550;
set(gcf, 'Position', pos);





