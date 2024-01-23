clear all; close all;

C1 = [204 247 255]./255 ; C11 = [80 180 255]./255 ; 
C3 = [179 179 255]./255; C33 = [25 25 255]./255;
C2 = [255 221 200]./255; C22 = [255 140 0]./255;
C4 = [255 200 200]./255; C44 = [200 0 0]./255; 
ms = 5; ls = 0.5; lw2 = 2; lw = 1.5;
xmin = 900; xmax = 2000; bin = 900;
pos(3) = 350; pos(4) = 1000;

% initial set up
load('sample9.mat');
dat5(3,2) = 5; dat5(5,2) = -4.5; dat5(4,2) = -3.5; 


x = [-500:1:4000];
TL = 1800 ; penalty = 0;
for i = 1:length(x)
    if x(i) < 0
        gainFunczero(i) = 0;
    elseif x(i) >= 0 && x(i) <= TL
        gainFunczero(i) = x(i)/(TL)*100;
    elseif x(i) > TL
        gainFunczero(i) = penalty;
    end
end

TL = 1800 ; penalty = -500;
for i = 1:length(x)
    if x(i) < 0
        gainFunclarge(i) = 0;
    elseif x(i) >= 0 && x(i) <= TL
        gainFunclarge(i) = x(i)/(TL)*100;
    elseif x(i) > TL
        gainFunclarge(i) = penalty;
    end
end

%% Gaussian
for con = 1:4
    x = [-500:1:4000];
    if con == 1
        dat(:,:,con) = dat30;
        mu_hut = mean(dat30(:,1)); sigma_hut = std(dat30(:,1))*10;
        gainFunc = gainFunczero;
    elseif con == 2
        dat(:,:,con) = [dat5; nan(25,2)];        
        mu_hut = mean(dat5(:,1)); sigma_hut = std(dat5(:,1))*10;
        gainFunc = gainFunczero;
    elseif con == 3
        dat(:,:,con) = dat30;        
        mu_hut = mean(dat30(:,1)); sigma_hut = std(dat30(:,1))*10;
        gainFunc = gainFunclarge;
    elseif con == 4
        dat(:,:,con) = [dat5; nan(25,2)];
        mu_hut = mean(dat5(:,1)); sigma_hut = std(dat5(:,1))*10;
        gainFunc = gainFunclarge;
    end
    
    for mu = 1:1:2500
        gaussian = exp(-(x - mu).^2/(2*sigma_hut^2)) / sqrt(2*pi*sigma_hut^2) ;
        gaussian = gaussian / sum(gaussian);
        EG_g(mu,con) = sum(gaussian.*gainFunc);        
    end
    
    x = 1:1:2500;
    [a d] = max(EG_g(:,con));
    OptAim_g(1,con) = d ;
    OptEG_g(1,con) = a;
    Sigmahut(1,con) = sigma_hut;    
end

for i = 1:4
    if i == 1
        C = C11; s = '-';
    elseif i == 2
        C = C22; s = '--';
    elseif i == 3
        C = C33; s = '-';
    elseif i == 4
        C = C44; s = '--';
    end
    
    subplot(10,1,i*2-1);
    plot(dat(:,1,i)*10+OptAim_g(1,i),dat(:,2,i)*10, 'ko','MarkerFaceColor','w','MarkerSize',ms,'linewidth',ls); 
    xlim([xmin xmax]); 
    yticks(-1000:2000:1000); xticks(-1000:4000:3000);
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1);
    
    subplot(10,1,i*2);
    gaussian = exp(-(x - OptAim_g(1,i)).^2/(2*Sigmahut(1,i)^2)) / sqrt(2*pi*Sigmahut(1,i)^2) ;
    gaussian = gaussian / sum(gaussian); 
    plot(x,gaussian,'-', 'Color',C, 'linewidth',lw); hold on 
    lineplot(OptAim_g(1,i), 'v','k-', 'linewidth',1);
    lineplot(1800, 'v','g-', 'linewidth',1);
    
    xlim([xmin xmax]); ylim([0 0.004]); 
    yticks(0:0.002:0.004); xticks(xmin:bin:xmax);
    xticklabels({'','',''}); yticklabels({'','',''});    
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize',1, 'linewidth', 1,'TickLength',[0.016 0]);

    subplot(10,1,9:10);
    plot(x, EG_g(:,i),s,'Color',C, 'linewidth',lw2); hold on
    ylim([-100 100]); xlim([xmin xmax]);
    xticks(xmin:bin:xmax);
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.016 0]);

end

set(gcf, 'Position', pos);

%% Laplacian
for con = 1:4
    x = [-500:1:4000];
    if con == 1
        dat(:,:,con) = dat30;
        median_hut = median(dat30(:,1));
        b_hut = mean(abs(dat30(:,1) - median_hut))*10;        
        gainFunc = gainFunczero;
    elseif con == 2
        dat(:,:,con) = [dat5; nan(25,2)];        
        median_hut = median(dat5(:,1));
        b_hut = mean(abs(dat5(:,1) - median_hut))*10;        
        gainFunc = gainFunczero;
    elseif con == 3
        dat(:,:,con) = dat30;        
        median_hut = median(dat30(:,1));
        b_hut = mean(abs(dat30(:,1) - median_hut))*10;        
        gainFunc = gainFunclarge;
    elseif con == 4
        dat(:,:,con) = [dat5; nan(25,2)];
        median_hut = median(dat5(:,1));
        b_hut = mean(abs(dat5(:,1) - median_hut))*10;        
        gainFunc = gainFunclarge;
    end
    
    eg_l = [];    
    for mu = 1:1:2500
        laplacian = 1/(2*b_hut) * exp(-(abs(x-mu)/b_hut)) ;
        laplacian = laplacian / sum(laplacian);
        eg_l = [eg_l; sum(laplacian.*gainFunc), mu];        
    end    
    EG_l(:,:,con) = eg_l;    
    x = 1:1:2500;
    [a d] = max(EG_l(:,1,con));    
    OptAim_l(1,con) = EG_l(d,2,con) ;
    OptEG_l(1,con) = a;
    Bhut(1,con) = b_hut;    
end

figure
for i = 1:4
    if i == 1
        C = C11; s = '-';
    elseif i == 2
        C = C22; s = '--';
    elseif i == 3
        C = C33; s = '-';
    elseif i == 4
        C = C44; s = '--';
    end
    
    subplot(10,1,i*2-1);
    plot(dat(:,1,i)*10+OptAim_l(1,i),dat(:,2,i)*10, 'ko','MarkerFaceColor','w','MarkerSize',ms,'linewidth',ls); 
    xlim([xmin xmax]); 
    yticks(-1000:2000:1000); xticks(-1000:4000:3000);
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1);
        
    subplot(10,1,i*2);
    laplacian = 1/(2*Bhut(1,i)) * exp(-(abs(x-OptAim_l(1,i))/Bhut(1,i))) ;
    laplacian = laplacian / sum(laplacian);    
    plot(x,laplacian,'-', 'Color',C, 'linewidth',lw); hold on 
    xlim([xmin xmax]); ylim([0 0.006]);    
    yticks(0:0.003:0.006); xticks(xmin:bin:xmax);
    xticklabels({'','',''}); yticklabels({'','',''});        
    lineplot(OptAim_l(1,i), 'v','k-', 'linewidth',1);
    lineplot(1800, 'v','g-', 'linewidth',1);    
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize',0.01, 'linewidth', 1,'TickLength',[0.016 0]);    

    subplot(10,1,9:10);
    plot(EG_l(:,2,i), EG_l(:,1,i),s,'Color',C, 'linewidth',lw2); hold on
    ylim([-100 100]); xlim([xmin xmax]);
    xticks(xmin:bin:xmax);
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.016 0]);
end
set(gcf, 'Position', pos);


%% Unifrom
for con = 1:4
    x = [-500:1:4000];
    if con == 1
        dat(:,:,con) = dat30;
        L_hut = min(dat30(:,1))*10;
        Smax = max(dat30(:,1))*10;        
        gainFunc = gainFunczero;
    elseif con == 2
        dat(:,:,con) = [dat5; nan(25,2)];        
        L_hut = min(dat5(:,1))*10;
        Smax = max(dat5(:,1))*10;        
        gainFunc = gainFunczero;
    elseif con == 3
        dat(:,:,con) = dat30;        
        L_hut = min(dat30(:,1))*10;
        Smax = max(dat30(:,1))*10;        
        gainFunc = gainFunclarge;
    elseif con == 4
        dat(:,:,con) = [dat5; nan(25,2)];
        L_hut = min(dat5(:,1))*10;
        Smax = max(dat5(:,1))*10;        
        gainFunc = gainFunclarge;
    end
    
    eg_u = [];
    for mu = 1:1:2500
        Mu = (L_hut + Smax)/2 + mu;
        uniform = zeros(1,length(x)) ;
        uniform(round(L_hut+mu)+abs(min(x)):round(Smax+mu)+abs(min(x))) = 1/(L_hut - Smax);
        uniform = uniform/sum(uniform) ; 
        eg_u = [eg_u; sum(uniform.*gainFunc), Mu];
    end
    EG_u(:,:,con) = eg_u;
    x = 1:1:2500;
    [a d] = max(EG_u(:,1,con));
    OptShift(1,con) = x(d) ;
    OptAim_u(1,con) = EG_u(d,2,con) ;
    OptEG_u(1,con) = a;
    Lhut(1,con) = L_hut; 
    Mhut(1,con) = Smax; 
end

figure
for i = 1:4
    if i == 1
        C = C11; s = '-';
    elseif i == 2
        C = C22; s = '--';
    elseif i == 3
        C = C33; s = '-';
    elseif i == 4
        C = C44; s = '--';
    end
    
    subplot(10,1,i*2-1);
    plot(dat(:,1,i)*10+OptShift(1,i),dat(:,2,i)*10, 'ko','MarkerFaceColor','w','MarkerSize',ms,'linewidth',ls); 
    xlim([xmin xmax]); 
    yticks(-1000:2000:1000); xticks(-1000:4000:3000);
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1);
    
    subplot(10,1,i*2);
    uniform = zeros(1,length(x)) ;
    uniform(round(OptShift(1,i)+Lhut(1,i))+abs(min(x)):round(OptShift(1,i)+Mhut(1,i))+abs(min(x))) = 1/(Lhut(1,i) - Mhut(1,i));
    uniform = uniform/sum(uniform) ;
        
    plot(x,uniform,'-', 'Color',C, 'linewidth',lw); hold on 
    xlim([xmin xmax]); ylim([0 0.003]);
    yticks(0:0.0015:0.003); xticks(xmin:bin:xmax);
    xticklabels({'','',''}); yticklabels({'','',''});        
    lineplot(OptAim_u(1,i), 'v','k-', 'linewidth',1);
    lineplot(1800, 'v','g-', 'linewidth',1);
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize',0.01, 'linewidth', 1,'TickLength',[0.016 0]);    

    subplot(10,1,9:10);
    plot(EG_u(:,2,i), EG_u(:,1,i),s,'Color',C, 'linewidth',lw2); hold on
    ylim([-100 100]); xlim([xmin xmax]);
    xticks(xmin:bin:xmax);
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.016 0]);
end

set(gcf, 'Position', pos);


%% Smax 
figure
for con = 1:4
    if con == 1
        dat(:,:,con) = dat30;
        Smax(1,con) = max(dat30(:,1))*10;        
        bias(1,con) = 50; 
    elseif con == 2
        dat(:,:,con) = [dat5; nan(25,2)];  
        Smax(1,con) = max(dat5(:,1))*10;        
        bias(1,con) = 50;         
    elseif con == 3
        dat(:,:,con) = dat30;        
        Smax(1,con) = max(dat30(:,1))*10;        
        bias(1,con) = 150;         
    elseif con == 4
        dat(:,:,con) = [dat5; nan(25,2)];
        Smax(1,con) = max(dat5(:,1))*10;        
        bias(1,con) = 150;         
    end
    
    yhut1(1,con) = (1800 - bias(1,con)) - Smax(:,con) ;
    
    subplot(8,1,con);
    plot(dat(:,1,con)*10+yhut1(1,con),dat(:,2,con)*10, 'ko','MarkerFaceColor','w','MarkerSize',ms,'linewidth',ls); hold on
    lineplot(1800, 'v','g-', 'linewidth',1);    
    xlim([xmin xmax]); 
    yticks(-1000:2000:1000); xticks(-1000:4000:3000);
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1);
end
set(gcf, 'Position', pos);

%% S50
figure
for con = 1:4
    if con == 1
        dat(:,:,con) = dat30;
        sdat = sort(dat30(:,1));         
        S50(1,con) = sdat(15)*10;        
        bias(1,con) = 300; 
    elseif con == 2
        dat(:,:,con) = [dat5; nan(25,2)];  
        sdat = sort(dat5(:,1));         
        S50(1,con) = sdat(3)*10;        
        bias(1,con) = 300;         
    elseif con == 3
        dat(:,:,con) = dat30;        
        sdat = sort(dat30(:,1));         
        S50(1,con) = sdat(15)*10;        
        bias(1,con) = 400;         
    elseif con == 4
        dat(:,:,con) = [dat5; nan(25,2)];
        sdat = sort(dat5(:,1));         
        S50(1,con) = sdat(3)*10;        
        bias(1,con) = 400;         
    end
    
    yhut2(1,con) = (1800 - bias(1,con)) - S50(:,con) ;
    
    subplot(8,1,con);
    plot(dat(:,1,con)*10+yhut2(1,con),dat(:,2,con)*10, 'ko','MarkerFaceColor','w','MarkerSize',ms,'linewidth',ls); hold on 
    lineplot(1800, 'v','g-', 'linewidth',1);    
    xlim([xmin xmax]); 
    yticks(-1000:2000:1000); xticks(-1000:4000:3000);
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1);
end
set(gcf, 'Position', pos);


%% Smin
figure
for con = 1:4
    if con == 1
        dat(:,:,con) = dat30;
        Smin(1,con) = min(dat30(:,1))*10;        
        bias(1,con) = 520; 
    elseif con == 2
        dat(:,:,con) = [dat5; nan(25,2)];  
        Smin(1,con) = min(dat5(:,1))*10;        
        bias(1,con) = 520;         
    elseif con == 3
        dat(:,:,con) = dat30;        
        Smin(1,con) = min(dat30(:,1))*10;        
        bias(1,con) = 650;         
    elseif con == 4
        dat(:,:,con) = [dat5; nan(25,2)];
        Smin(1,con) = min(dat5(:,1))*10;        
        bias(1,con) = 650;         
    end
    
    yhut3(1,con) = (1800 - bias(1,con)) - Smin(:,con) ;
    
    subplot(8,1,con);
    plot(dat(:,1,con)*10+yhut3(1,con),dat(:,2,con)*10, 'ko','MarkerFaceColor','w','MarkerSize',ms,'linewidth',ls); hold on 
    lineplot(1800, 'v','g-', 'linewidth',1);    
    xlim([xmin xmax]); 
    yticks(-1000:2000:1000); xticks(-1000:4000:3000);
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1);
end
set(gcf, 'Position', pos);

