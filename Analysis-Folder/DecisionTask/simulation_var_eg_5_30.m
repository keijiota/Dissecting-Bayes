clear all; close all; 

load var2optaim2
var2optaim = var2optaim2 ;

x = [-100:0.1:300];
TL = 180 ; penalty = 0;
for i = 1:length(x)
    if x(i) < 0
        gainFunc(i) = 0;
    elseif x(i) >= 0 && x(i) <= TL
        gainFunc(i) = x(i)/(TL)*100;
    elseif x(i) > TL
        gainFunc(i) = penalty;
    end
end

N1 = 5; N2 = 30; Nrep = 10000;
Mu = [0; 0]; 
sigma_y = 10; %10~20 mm 
Cov = [sigma_y^2/10 0; 0 sigma_y^2];
for irep = 1:Nrep
sample1 = ndRandn(Mu,Cov,N1); sample1 = sample1';
sample2 = ndRandn(Mu,Cov,N2); sample2 = sample2';

% sufficient statistics
samplemean1 = mean(sample1); samplecov1  = cov(sample1);
samplemean2 = mean(sample1); samplecov2  = cov(sample2);
samplevar_y(irep,1:2) = [samplecov1(2,2), samplecov2(2,2)]; 

% from sufficient statistics to ideal observer 
% EG_g = [];
% for mu = 100:0.1:180
%     gaussian = exp(-(x- mu).^2./(2*samplevar_y(irep,1:2)')) ./ sqrt(2*pi*samplevar_y(irep,1:2)') ; 
%     gaussian = gaussian ./ sum(gaussian,2);
%     EG_g = [EG_g; sum(gaussian.*gainFunc,2)', mu];    
% end
% 
% [opteg(irep,1) tmp] = max(EG_g(:,1));
% optaim(irep,1) = EG_g(tmp,3) ;
% [opteg(irep,2) tmp] = max(EG_g(:,2));
% optaim(irep,2) = EG_g(tmp,3) ;
end

% optaim = interp1(linspace(0,4000,4000), linspace(180,0,4000),samplevar_y); 
optaim = interp1(var2optaim(:,3), var2optaim(:,1), samplevar_y);


figure
bw = 50; d = 0.4; ymin = 90;
subplot(2,3,1)
histogram(samplevar_y(:,1), 'FaceAlpha', d, 'FaceColor','b','Normalization', 'probability', 'BinWidth', bw); hold on 
histogram(samplevar_y(:,2), 'FaceAlpha', d, 'FaceColor','r','Normalization', 'probability', 'BinWidth', bw);
lineplot(mean(samplevar_y(:,1)), 'v', 'b', 'linewidth', 2);
lineplot(mean(samplevar_y(:,2)), 'v', 'r--', 'linewidth', 2); xlim([0 1000]); 
xlabel('Sample variance'); ylabel('Normalized histogram'); 
title('Population variance_y is fixed at 225 mm');legend({'N = 5','N = 30'}); 

bw = 2.5;
subplot(2,3,2)
histogram(optaim(:,1), 'FaceAlpha', d, 'FaceColor','b','Normalization', 'probability', 'BinWidth', bw); hold on 
histogram(optaim(:,2), 'FaceAlpha', d, 'FaceColor','r','Normalization', 'probability', 'BinWidth', bw);
lineplot(nanmean(optaim(:,1)), 'v', 'b', 'linewidth', 2);
lineplot(nanmean(optaim(:,2)), 'v', 'r--', 'linewidth', 2);
xlabel('Ideal set point'); ylabel('Normalized histogram'); xlim([ymin 180]); 
title(strcat('Blue solid line=', num2str(nanmean(optaim(:,1))), ', Red dashed line', num2str(nanmean(optaim(:,2)))));

subplot(2,3,3)
plot(samplevar_y, optaim,'.');
xlim([ymin 180]); xlim([0 1000]); 
xlabel('Sample variance'); ylabel('Ideal set point');
title('Ideal set point as a function of sample variance');legend({'N = 5','N = 30'}); 

%% random covariance

sigma_y = rand(Nrep,1)*10+10; %10~20 mm 

for irep = 1:Nrep
sigma_y = rand(1,1)*10+10 ;   
Cov = [sigma_y^2/10 0; 0 sigma_y^2];    
sample11 = ndRandn(Mu,Cov,N1); sample11 = sample11';
sample22 = ndRandn(Mu,Cov,N2); sample22 = sample22';

% sufficient statistics
samplemean11 = mean(sample11); samplecov11  = cov(sample11);
samplemean22 = mean(sample11); samplecov22  = cov(sample22);
samplevar_y2(irep,1:2) = [samplecov11(2,2), samplecov22(2,2)]; 
end

% from sufficient statistics to ideal observer 
optaim2(:,1) = interp1(var2optaim(:,3), var2optaim(:,1), samplevar_y2(:,1));
optaim2(:,2) = interp1(var2optaim(:,3), var2optaim(:,1), samplevar_y2(:,2));

bw = 50; d = 0.4;
subplot(2,3,4)
histogram(samplevar_y2(:,1), 'FaceAlpha', d, 'FaceColor','b','Normalization', 'probability', 'BinWidth', bw); hold on 
histogram(samplevar_y2(:,2), 'FaceAlpha', d, 'FaceColor','r','Normalization', 'probability', 'BinWidth', bw);
lineplot(mean(samplevar_y2(:,1)), 'v', 'b', 'linewidth', 2);
lineplot(mean(samplevar_y2(:,2)), 'v', 'r--', 'linewidth', 2);
xlabel('Sample variance'); ylabel('Normalized histogram'); xlim([0 1500]); 
title('Population variance_y is randomly sampled from 100 to 400 mm'); legend({'N = 5','N = 30'}); 

bw = 2.5;
subplot(2,3,5)
histogram(optaim2(:,1), 'FaceAlpha', d, 'FaceColor','b','Normalization', 'probability', 'BinWidth', bw); hold on 
histogram(optaim2(:,2), 'FaceAlpha', d, 'FaceColor','r','Normalization', 'probability', 'BinWidth', bw);
lineplot(nanmean(optaim2(:,1)), 'v', 'b', 'linewidth', 2);
lineplot(nanmean(optaim2(:,2)), 'v', 'r--', 'linewidth', 2);
xlabel('Ideal set point'); ylabel('Normalized histogram'); xlim([ymin 180]); 
title(strcat('Blue solid line=', num2str(nanmean(optaim2(:,1))), ', Red dashed line', num2str(nanmean(optaim2(:,2)))));

subplot(2,3,6)
plot(samplevar_y2,optaim2,'.');
xlabel('Sample variance'); ylabel('Ideal set point');
title('Ideal set point as a function of sample variance');legend({'N = 5','N = 30'}); 
xlim([ymin 180]); xlim([0 1500]); 




