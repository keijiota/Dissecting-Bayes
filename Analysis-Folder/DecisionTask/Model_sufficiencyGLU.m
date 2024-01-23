clear all; close all;

% internal model from sample
% MLE of mu in a Gaussian -> Exp(Data)
% MLE of sigma in a Gaussian -> sum(D - mu)^2 / (n-1)
% MLE of mu in a Laplacian -> sample median
% MLE of b in a Laplacian -> sum(abs(D - mu))/ n, mean asolute deviation from median 
% MLE of upper and lower bound in a uniform -> min(D), max(D)

x = -5:0.01:5;
sigma = 1; N = 5; penalty = -500;
sample = randn(N,1) * sigma ; 

mu_hut = mean(sample);
sigma_hut = std(sample);
median_hut = median(sample);
b_hut = mean(abs(sample - median_hut));
L_hut = min(sample);
M_hut = max(sample);

gaussian = exp(-(x - mu_hut).^2/(2*sigma_hut^2)) / sqrt(2*pi*sigma_hut^2) ;
gaussian = gaussian / sum(gaussian);

laplacian = 1/(2*b_hut) * exp(-(abs(x-median_hut)/b_hut)) ;
laplacian = laplacian / sum(laplacian);

uniform = zeros(1,length(x)) ;
uniform(round(length(x)/2+L_hut*100):round(length(x)/2+M_hut*100)) = 1/(L_hut - M_hut);
uniform = uniform/sum(uniform) ;

figure
subplot(2,1,1)
plot(sample, ones(N,1),'ko'); 
xlim([min(x) max(x)]);

subplot(2,1,2)
plot(x,gaussian); hold on 
plot(x,laplacian); 
plot(x,uniform);


x = [-500:1:4000];
TL = 1800 ;
for i = 1:length(x)
    if x(i) < 0
       gainFunc(i) = 0;        
    elseif x(i) >= 0 && x(i) <= TL 
       gainFunc(i) = x(i)/(TL)*100;
    elseif x(i) > TL
       gainFunc(i) = penalty;
    end    
end

sigma = 150;
sample = randn(N,1) * sigma ; 

mu_hut = mean(sample);
sigma_hut = std(sample);
median_hut = median(sample);
b_hut = mean(abs(sample - median_hut));
L_hut = min(sample);
M_hut = max(sample);

EG_g = []; EG_l = []; EG_u = [];
for mu = 1:1:2500
    gaussian = exp(-(x - mu).^2/(2*sigma_hut^2)) / sqrt(2*pi*sigma_hut^2) ;
    gaussian = gaussian / sum(gaussian);
    EG_g = [EG_g; sum(gaussian.*gainFunc), mu];
    
    laplacian = 1/(2*b_hut) * exp(-(abs(x-mu)/b_hut)) ;
    laplacian = laplacian / sum(laplacian);
    EG_l = [EG_l; sum(laplacian.*gainFunc), mu];

    Mu = (L_hut + M_hut)/2 + mu;    
    uniform = zeros(1,length(x)) ;
    uniform(round(L_hut+mu)+abs(min(x)):round(M_hut+mu)+abs(min(x))) = 1/(L_hut - M_hut);
    uniform = uniform/sum(uniform) ;
    EG_u = [EG_u; sum(uniform.*gainFunc), Mu];
end

x = 1:1:2500;
[OptEG_g d] = max(EG_g(:,1));
OptAim_g = EG_g(d,2) ;
gaussian = exp(-(x - OptAim_g).^2/(2*sigma_hut^2)) / sqrt(2*pi*sigma_hut^2) ;
gaussian = gaussian / sum(gaussian);

[OptEG_l d] = max(EG_l(:,1));
OptAim_l = EG_l(d,2) ;
laplacian = 1/(2*b_hut) * exp(-(abs(x-OptAim_l)/b_hut)) ;
laplacian = laplacian / sum(laplacian);

[OptEG_u d] = max(EG_u(:,1));
OptShift = x(d) ;
OptAim_u = EG_u(d,2) ;
uniform = zeros(1,length(x)) ;
uniform(round(OptShift+L_hut)+abs(min(x)):round(OptShift+M_hut)+abs(min(x))) = 1/(L_hut - M_hut);
uniform = uniform/sum(uniform) ;

figure
subplot(3,3,1)
plot(x,gainFunc(501:3000),'r-');
subplot(3,3,4)
plot(x,gaussian,'k-');
subplot(3,3,7)
plot(EG_g(:,2),EG_g(:,1),'k-'); hold on 
lineplot(OptAim_g, 'v','k-'); 
xlim([min(x) max(x)]);

subplot(3,3,2)
plot(x,gainFunc(501:3000),'r-');
subplot(3,3,5)
plot(x,laplacian,'k-');
subplot(3,3,8)
plot(EG_l(:,2),EG_l(:,1),'k-'); hold on 
lineplot(OptAim_l, 'v','k-'); 
xlim([min(x) max(x)]);

subplot(3,3,3)
plot(x,gainFunc(501:3000),'r-');
subplot(3,3,6)
plot(x,uniform,'k-');
lineplot(OptAim_u, 'v','k-'); xlim([min(x) max(x)]);
subplot(3,3,9)
plot(EG_u(:,2),EG_u(:,1),'k-'); hold on 
lineplot(OptAim_u, 'v','k-'); 
xlim([min(x) max(x)]);

