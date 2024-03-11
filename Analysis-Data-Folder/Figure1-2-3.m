clear; close all;

x = -5:.01:5; y = -5:.01:5;
sdx = 0.75; sdy = 1.5; r = 0;
a = 550; b = 550; ms = 11;
c = jet; c = c(1:200,:);
% c = parula; c = c(1:end,:);

%% Fig.1
%------Fig.1A----------------
[X,Y] = meshgrid(x,y);
mux = -1; muy = -1;
fxy = 1/(2*pi*sdx*sdy*sqrt(1-(r^2)))...
    * exp(-1/(1-(r^2)) * (((X-mux).^2 /2/sdx^2) + ((Y-muy).^2 /2/sdy^2) - (r*(X-mux).*(Y-muy)/sdx/sdy))) * diff(x(1:2)) * diff(y(1:2));

sum(sum(fxy))

figure;
mesh(X,Y,fxy); axis('equal'); colormap(c); 
xlim([-2.7 2.7]); ylim([-4.5 4.5]); xticks(-10:20:10); yticks(-10:20:10);
view(2);
pos(3) = a; pos(4) = b;
set(gcf, 'Position', pos); 
caxis([0 0.000014]);

%------Fig.1B----------------
mux = 0; muy = 0;
fxy = 1/(2*pi*sdx*sdy*sqrt(1-(r^2)))...
    * exp(-1/(1-(r^2)) * (((X-mux).^2 /2/sdx^2) + ((Y-muy).^2 /2/sdy^2) - (r*(X-mux).*(Y-muy)/sdx/sdy))) * diff(x(1:2)) * diff(y(1:2));

figure;
mesh(X,Y,fxy); axis('equal'); colormap(c); 
xlim([-2.7 2.7]); ylim([-4.5 4.5]); xticks(-10:20:10); yticks(-10:20:10);
view(2);
pos(3) = a; pos(4) = b;
set(gcf, 'Position', pos);
caxis([0 0.000014]);
 
%------Fig.1C----------------
mux = 0; muy = 0;
fxy = 1/(2*pi*sdx*sdy*sqrt(1-(r^2)))...
    * exp(-1/(1-(r^2)) * (((X-mux).^2 /2/sdx^2) + ((Y-muy).^2 /2/sdy^2) - (r*(X-mux).*(Y-muy)/sdx/sdy))) * diff(x(1:2)) * diff(y(1:2));

figure
mesh(X,Y,fxy); axis('equal'); colormap(c); 
xlim([-2.7 2.7]); ylim([-4.5 4.5]); xticks(-10:20:10); yticks(-10:20:10);
view(2);
pos(3) = a; pos(4) = b;
set(gcf, 'Position', pos);
caxis([0 0.000014]);
 
%% Fig.2
%------Fig.2A----------------
mux = 0; muy = 0;
fxy = 1/(2*pi*sdx*sdy*sqrt(1-(r^2)))...
    * exp(-1/(1-(r^2)) * (((X-mux).^2 /2/sdx^2) + ((Y-muy).^2 /2/sdy^2) - (r*(X-mux).*(Y-muy)/sdx/sdy))) * diff(x(1:2)) * diff(y(1:2));

figure
mesh(X,Y,fxy); axis('equal'); colormap(c); 
xlim([-2.7 2.7]); ylim([-4.5 4.5]); xticks(-10:20:10); yticks(-10:20:10);
view(2);
pos(3) = a; pos(4) = b;
set(gcf, 'Position', pos);
caxis([0 0.000014]);
 
%------Fig.2B----------------
n = 30 ;
meanmat = [mux, muy]';
covmat = [sdx^2, r; r, sdy^2];

runsampling = 0;
if runsampling == 1
    mux_t = 0.2; muy_t = 0.2;
    sdx_t = 1.5; sdy_t = 3;
    raw_t = 0.15;
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
end

load('s30_4');
s = S(:,:,10);
musx = mean(s(1,:)); musy = mean(s(2,:));
covsxy = cov(s') ;
sdsx = covsxy(1)^(1/2); sdsy = covsxy(4)^(1/2); rawxy = covsxy(2);

figure
plot(s(1,:), s(2,:),'ko','MarkerFaceColor','w', 'MarkerSize',ms);
axis('equal');
xlim([-2.7 2.7]); ylim([-4.5 4.5]);
xticks(-10:20:10); yticks(-10:20:10);
pos(3) = a; pos(4) = b;
set(gcf, 'Position', pos);

%% Fig.3
%------Fig.3A&B----------------
mux = 0; muy = 0;
fxy = 1/(2*pi*sdx*sdy*sqrt(1-(r^2)))...
      * exp(-1/(1-(r^2)) * (((X-mux).^2 /2/sdx^2) + ((Y-muy).^2 /2/sdy^2) - (r*(X-mux).*(Y-muy)/sdx/sdy))) * diff(x(1:2)) * diff(y(1:2));

figure
mesh(X,Y,fxy); axis('equal'); colormap(c); 
xlim([-2.7 2.7]); ylim([-4.5 4.5]); xticks(-10:20:10); yticks(-10:20:10);
view(2);
pos(3) = a; pos(4) = b;
set(gcf, 'Position', pos);  
% figure
% c = colorbar; %caxis([0 0.000014]);

figure
contour(X,Y,fxy,7,'-','linewidth',0.75); hold on
axis('equal');
xlim([-2.7 2.7]); ylim([-4.5 4.5]); xticks(-10:20:10); yticks(-10:20:10);
view(2);  
caxis([0 0.000014]);

load('s8_1'); 
[sx tmp] = sort(s(1,:)); sy = s(2,tmp);
s = [sx;sy];
s(1,4) = s(1,4)-0.5;
s(1,5) = s(1,5)+0.25;
s(2,8) = s(2,8)-0.5; s(1,8) = s(1,8)-0.2;
s2 = [s(1,:); -s(2,:)];

mus = mean(s,2)
covsxy = cov(s'); [covsxy(1,1)^(1/2) covsxy(1,2) covsxy(2,2)^(1/2)]

plot(s(1,:), s(2,:),'ko','MarkerFaceColor','w', 'MarkerSize',ms); hold on
axis('equal');
xlim([-2.7 2.7]); ylim([-4.5 4.5]); xticks(-10:20:10); yticks(-10:20:10);
pos(3) = a; pos(4) = b;
set(gcf, 'Position', pos);

 








