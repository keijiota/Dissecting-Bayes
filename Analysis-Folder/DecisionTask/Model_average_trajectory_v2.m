clear all; close all ;

load('dat_decision');

x = [-1000:1:3000];
TL = 1800 ;
for i = 1:length(x)
    if x(i) < 0
        gainFunc(i) = 0;
    elseif x(i) >= 0 && x(i) <= TL
        gainFunc(i) = x(i)/(TL)*100;
    elseif x(i) > TL
        gainFunc(i) = 0;
    end
end
for i = 1:length(x)
    if x(i) < 0
        gainFunc2(i) = 0;
    elseif x(i) >= 0 && x(i) <= TL
        gainFunc2(i) = x(i)/(TL)*100;
    elseif x(i) > TL
        gainFunc2(i) = -500;
    end
end

% making a matrix we need
for condition = 1:4
    if condition == 1 || condition == 3
        gainfunc = gainFunc;
    else
        gainfunc = gainFunc2;
    end
    
    endpointorg = EndpointorgY(:,:,condition);
    samplesY = SamplesY(:,:,condition);
    orgnoise = OrgNoiseY(:,:,condition);
    samplenoise = SampleNoiseY(:,:,condition);
    optAim = OptAim(:,:,condition);
    ssy = []; sy = []; ssamples = [];
    
    [N T] = size(endpointorg);
    for subi = 1:N
        for itrial = 1:T
            % each samples sorted in a order
            sy = samplesY(subi, itrial*31-30:itrial*31-1) ; %X1-X30
            endpoint = endpointorg(subi, :);

            if condition == 3 || condition == 4
                sy = sy(1:5);
            end
            ssy = sort(sy);
            if condition == 1 || condition == 2
                pctl1 = round(30*0.02); pctl2 = round(30*0.25); pctl3 = round(30*0.5);
                pctl4 = round(30*0.75); pctl5 = round(30*1);
                ssy = [ssy(pctl1), ssy(pctl2),ssy(pctl3),ssy(pctl4),ssy(pctl5)];
            end
            ssamples(itrial, :, subi) = ssy;
            
            % internal model from sample
            % MLE of mu in a Gaussian -> Exp(Data)
            % MLE of sigma in a Gaussian -> sum(D - mu)^2 / (n-1)
            % MLE of mu in a Laplacian -> sample median
            % MLE of b in a Laplacian -> sum(abs(D - mu))/ n, mean asolute deviation from median
            % MLE of upper and lower bound in a uniform -> min(D), max(D)
            mu_hut = mean(sy);
            Sigma_hut(subi,itrial) = std(sy);
            median_hut = median(sy);
            B_hut(subi,itrial) = mean(abs(sy - median_hut));
            L_hut(subi,itrial) = min(sy - endpoint(itrial));
            M_hut(subi,itrial) = max(sy - endpoint(itrial));
        end
    end
    
    for subi = 1:N
        onoise = orgnoise(subi,:)';
        snoise = samplenoise(subi,:)';
        endpoint = endpointorg(subi, :)';
        optaim = optAim(subi,:)';
        sample = ssamples(:,:,subi);
        S = (sample-endpoint)./onoise;
        orgsample = sample - endpoint ;
        
        normSmin(:, condition, subi) = S(:,1) ;
        normSmax(:, condition, subi) = S(:,5) ;
        Smin(:, condition, subi) = orgsample(:,1) ;
        S50(:, condition, subi) = orgsample(:,3) ;        
        Smax(:, condition, subi) = orgsample(:,5) ;
        
        Endpoint(:, condition, subi) = endpoint ;
        Optaim(:, condition, subi) = optaim ;
        
        % argmax G(x)*P(x), P(x)~N(sigma_hut|sample), L(b_hut|sample), U(l_hut, m_hut|sample) 
        sigmahut = Sigma_hut(subi,:)';
        bhut = B_hut(subi,:)';
        lhut = L_hut(subi,:)';
        mhut = M_hut(subi,:)';
        
%         for itrial = 1:T
%             sigma_hut = sigmahut(itrial)*10; b_hut = bhut(itrial)*10;
%             l_hut = lhut(itrial)*10; m_hut = mhut(itrial)*10;            
%             EG_g = []; EG_l = []; EG_u = [];
%             for mu = 1:1:2000
%                 gaussian = exp(-(x - mu).^2/(2*sigma_hut^2)) / sqrt(2*pi*sigma_hut^2) ;
%                 gaussian = gaussian / sum(gaussian);
%                 EG_g = [EG_g; sum(gaussian.*gainfunc), mu];
%                 
%                 laplacian = 1/(2*b_hut) * exp(-(abs(x-mu)/b_hut)) ;
%                 laplacian = laplacian / sum(laplacian);
%                 EG_l = [EG_l; sum(laplacian.*gainfunc), mu];
%                 
%                 Mu = (l_hut + m_hut)/2 + mu;
%                 uniform = zeros(1,length(x)) ;
%                 uniform(round(l_hut+mu)+abs(min(x)):round(m_hut+mu)+abs(min(x))) = 1/(l_hut - m_hut);
%                 uniform = uniform/sum(uniform) ;
%                 EG_u = [EG_u; sum(uniform.*gainfunc), Mu];
%             end
%             
%             xx = 1:1:2000;
%             [OptEG_g d] = max(EG_g(:,1));
%             OptAim_g(itrial, condition, subi) = EG_g(d,2)/10 ;
%             
%             [OptEG_l d] = max(EG_l(:,1));
%             OptAim_l(itrial, condition, subi) = EG_l(d,2)/10 ;
%             
%             [OptEG_u d] = max(EG_u(:,1));
%             OptShift = xx(d);
%             OptAim_u(itrial, condition, subi) = EG_u(d,2)/10 ;
%         end        

    end % subject
end % condition

% save suffdata OptAim_g OptAim_l OptAim_u
load suffdata


tn = {'30 dots & 0','30 dots & -500','5 dots & 0','5 dots & -500'};

% model fitting
for subi = 1:N
    Optaim_g = OptAim_g(:,:,subi) ;
    Optaim_l = OptAim_l(:,:,subi) ;
    Optaim_u = OptAim_u(:,:,subi) ;       
    Subdata = Endpoint(:,:,subi) ;
    smax = Smax(:,:,subi) ;
    s50 = S50(:,:,subi) ;
    smin = Smin(:,:,subi) ;
    
    [estimates_smax, var_hut, log_likelihood] = ML_fit_SmaxSt(Subdata, smax);
    [estimates_s50, var_hut, log_likelihood] = ML_fit_S50St(Subdata, s50);
    [estimates_smin, var_hut, log_likelihood] = ML_fit_SminSt(Subdata, smin);    

    Bias1(subi,:) = estimates_smax;
    Bias2(subi,:) = estimates_s50;
    Bias3(subi,:) = estimates_smin;
    Yhut1 = []; Yhut2 = []; Yhut3 = [];    
    
%     figure
    for con = 1:4
        if con == 1 || con == 3
            bias1 = estimates_smax(1);
            bias2 = estimates_s50(1);            
            bias3 = estimates_smin(1);
        else
            bias1 = estimates_smax(2);
            bias2 = estimates_s50(2);                        
            bias3 = estimates_smin(2);
        end
        Yhut1(:,con) = (180 - bias1) - smax(:,con) ;
        Yhut2(:,con) = (180 - bias2) - s50(:,con) ;        
        Yhut3(:,con) = (180 - bias3) - smin(:,con) ;
        
%         subplot(2,2,con)
%         plot(Subdata(:,con), 'k-','linewidth',2); hold on
%         plot(Optaim_g(:,con), '-','linewidth',0.5,'Color',[1 0 0]); 
%         plot(Optaim_l(:,con), '-','linewidth',0.5, 'Color',[0 1 0]); 
%         plot(Optaim_u(:,con), '-','linewidth',0.5, 'Color',[0 0 1]);         
%         plot(Yhut1(:,con), '--','linewidth',1.5,'Color',[1 0 0]);
%         plot(Yhut2(:,con), '--','linewidth',1.5,'Color',[0 1 0]);        
%         plot(Yhut3(:,con), '--','linewidth',0.5, 'Color',[0 0 1]);         
%         title(tn(con));
%         if con == 4
%            legend({'Data','N(sigma)','L(b)','U(a,b)','Smax','S50','Smin'}); 
%         end
    end
    
    YHut1(:,:,subi) = Yhut1;
    YHut2(:,:,subi) = Yhut2;
    YHut3(:,:,subi) = Yhut3;    
    Estimates_smax(subi,:) = estimates_smax;
    Estimates_s50(subi,:) = estimates_s50;    
    Estimates_smin(subi,:) = estimates_smin;
    
    for con = 1:4
        subdata = Subdata(:,con);
        optaim_g = Optaim_g(:,con);
        optaim_l = Optaim_l(:,con);
        optaim_u = Optaim_u(:,con);
        yhut1 = Yhut1(:,con);
        yhut2 = Yhut2(:,con);
        yhut3 = Yhut3(:,con);
        [T M] = size(subdata);
        
        % argmax G(x)*P(x), P(x)~N(sigma_hut|sample)
        var_hut = mean((subdata - optaim_g).^2);
        log_likelihood = (T*M)/2 * (log(2*pi) + log(var_hut) + 1) ;
        K = 0; c = 1;
        AICc(subi, c, con) = -2*(-log_likelihood) + 2*K + 2*(K*(K+1))/(T*M-K-1) ;
        R2(subi, c, con) = 1 - (sum((optaim_g - subdata).^2) / sum((mean(subdata) - subdata).^2)) ;
        R2adj(subi, c, con) = 1 - ((1-R2(subi, c, con))*(T*M - 1) / (T*M - K - 1)) ;
        
        % argmax G(x)*P(x), P(x)~L(b_hut|sample)
        var_hut = mean((subdata - optaim_l).^2);
        log_likelihood = (T*M)/2 * (log(2*pi) + log(var_hut) + 1) ;
        K = 0; c = 2;
        AICc(subi, c, con) = -2*(-log_likelihood) + 2*K + 2*(K*(K+1))/(T*M-K-1) ;
        R2(subi, c, con) = 1 - (sum((optaim_l - subdata).^2) / sum((mean(subdata) - subdata).^2)) ;
        R2adj(subi, c, con) = 1 - ((1-R2(subi, c, con))*(T*M - 1) / (T*M - K - 1)) ;
        
        % argmax G(x)*P(x), P(x)~U(l_hut, m_hut|sample)
        var_hut = mean((subdata - optaim_u).^2);
        log_likelihood = (T*M)/2 * (log(2*pi) + log(var_hut) + 1) ;
        K = 0; c = 3;
        AICc(subi, c, con) = -2*(-log_likelihood) + 2*K + 2*(K*(K+1))/(T*M-K-1) ;
        R2(subi, c, con) = 1 - (sum((optaim_u - subdata).^2) / sum((mean(subdata) - subdata).^2)) ;
        R2adj(subi, c, con) = 1 - ((1-R2(subi, c, con))*(T*M - 1) / (T*M - K - 1)) ;
        
        % y = (PB - bias_gf) - Smax
        var_hut = mean((subdata - yhut1).^2);
        log_likelihood = (T*M)/2 * (log(2*pi) + log(var_hut) + 1) ;
        K = 1; c = 4;
        AICc(subi, c, con) = -2*(-log_likelihood) + 2*K + 2*(K*(K+1))/(T*M-K-1) ;
        R2(subi, c, con) = 1 - (sum((yhut1 - subdata).^2) / sum((mean(subdata) - subdata).^2)) ;
        R2adj(subi, c, con) = 1 - ((1-R2(subi, c, con))*(T*M - 1) / (T*M - K - 1)) ;
        
        % y = (PB - bias_gf) - S50
        var_hut = mean((subdata - yhut2).^2);
        log_likelihood = (T*M)/2 * (log(2*pi) + log(var_hut) + 1) ;
        K = 1; c = 5;
        AICc(subi, c, con) = -2*(-log_likelihood) + 2*K + 2*(K*(K+1))/(T*M-K-1) ;
        R2(subi, c, con) = 1 - (sum((yhut2 - subdata).^2) / sum((mean(subdata) - subdata).^2)) ;
        R2adj(subi, c, con) = 1 - ((1-R2(subi, c, con))*(T*M - 1) / (T*M - K - 1)) ;
        
        % y = (PB - bias_gf) - Smin
        var_hut = mean((subdata - yhut3).^2);
        log_likelihood = (T*M)/2 * (log(2*pi) + log(var_hut) + 1) ;
        K = 1; c = 6;
        AICc(subi, c, con) = -2*(-log_likelihood) + 2*K + 2*(K*(K+1))/(T*M-K-1) ;
        R2(subi, c, con) = 1 - (sum((yhut3 - subdata).^2) / sum((mean(subdata) - subdata).^2)) ;
        R2adj(subi, c, con) = 1 - ((1-R2(subi, c, con))*(T*M - 1) / (T*M - K - 1)) ;
    end
end

for subi = 1:N
    avSubdata(subi,:) = mean(Endpoint(:,:,subi));    
    avOptaim_g(subi,:) = mean(OptAim_g(:,:,subi));    
    avOptaim_l(subi,:) = mean(OptAim_l(:,:,subi));
    avOptaim_u(subi,:) = mean(OptAim_u(:,:,subi));
    avYHut1(subi,:) = mean(YHut1(:,:,subi));
    avYHut2(subi,:) = mean(YHut2(:,:,subi));
    avYHut3(subi,:) = mean(YHut3(:,:,subi));        
end

tmp = [1 3 2 4]; % [30,0], [5,0], [30,-500], [5,-500]
avSubdata = avSubdata(:,tmp);
avOptaim_g = avOptaim_g(:,tmp);
avOptaim_l = avOptaim_l(:,tmp);
avOptaim_u = avOptaim_u(:,tmp);
avYHut1 = avYHut1(:,tmp);
avYHut2 = avYHut2(:,tmp);
avYHut3 = avYHut3(:,tmp);
AICc = AICc(:,:,tmp) ; 
R2adj = R2adj(:,:,tmp) ; 

C1 = [204 247 255]./255 ; C11 = [80 180 255]./255 ; 
C3 = [179 179 255]./255; C33 = [25 25 255]./255;
C2 = [255 221 200]./255; C22 = [255 140 0]./255;
C4 = [255 200 200]./255; C44 = [200 0 0]./255; 
cs = 9; ms = 8; ms2 = 10;

figure
for iplot = [1:2, 5:7]
    if iplot == 1
        dat = avSubdata; yn = 'Endpoint [mm]';
    elseif iplot == 2
        dat = avOptaim_g; yn = 'Ehut by N(sigma)';
    elseif iplot == 3
        dat = avOptaim_l; yn = 'Ehut by L(b)';
    elseif iplot == 4
        dat = avOptaim_u; yn = 'Ehut by U(a,b)';
    elseif iplot == 5
        dat = avYHut1; yn = 'Ehut by Smax';
    elseif iplot == 6
        dat = avYHut2; yn = 'Ehut by S50';
    elseif iplot == 7
        dat = avYHut3; yn = 'Ehut by Smin';
    end
    
        for con = 1:4
            if con == 1
                C = C1; CC = C11; c = 0;
            elseif con == 2
                C = C2; CC = C22; c = 0.2;
            elseif con == 3
                C = C3; CC = C33; c = 0;
            elseif con == 4
                C = C4; CC = C44; c = 0.2;
            end
            subplot(3,4,iplot)
            plot(con, dat(:,con), 'o', 'MarkerEdgeColor', C, 'MarkerFaceColor', C, 'MarkerSize',ms); hold on
            plot(con, mean(dat(:,con)), 's', 'MarkerEdgeColor', CC, 'MarkerFaceColor', CC, 'MarkerSize',ms2); hold on
            errorbar(con, mean(dat(:,con)), std(dat(:,con)), 'o', 'color', CC, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
        end
        plot(1:2, [mean(dat(:,1)), mean(dat(:,2))], 'k--', 'linewidth', 1);
        plot(3:4, [mean(dat(:,3)), mean(dat(:,4))], 'k--', 'linewidth', 1);
        ylim([110 180]); xlim([0.25 4.75]); xticks(1:1:8); yticks(100:20:180);
        xticklabels({'30, 0', '5, 0', '30, -500', '5, -500'});
        ylabel(yn, 'FontName', 'Arial', 'FontSize', 10);
        set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);
end

subplot(3,4,10:11)
nmodel = 4;
for con = 1:4
    if con == 1
        C = C11; c = 0;
    elseif con == 2
        C = C22; c = 0.2;
    elseif con == 3
        C = C33; c = 0.1;
    elseif con == 4
        C = C44; c = 0.3;
    end    
    plot(1+c:nmodel+c, mean(AICc(:,[1,4:6],con)), '-s', 'Color', C, 'MarkerEdgeColor', C, 'MarkerFaceColor', C, 'linewidth',1,'MarkerSize',ms2); hold on
%     errorbar(1+c:nmodel+c, mean(AICc(:,:,con)), std(AICc(:,:,con)), 'o', 'color', C, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
    ylabel('AICc');
    xticks(1+0.15:1:nmodel+0.15); xlim([0.25 nmodel+0.75]);
    ylim([340 460]); yticks(340:40:460);
    set(gca,'XTickLabel',{'N(sigma)','L(b)','U(a,b)','Smax','S50','Smin'});
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.016 0]);
    legend({'30 dots, 0 penalty', '','5 dots, 0 penalty', '','30 dots, -500 penalty', '','5 dots, -500 penalty'});    
end
r = 1500/1600;
% pos(3) = 900; pos(4) = r*900;
pos(3) = 1500; pos(4) = 1600;
set(gcf, 'Position', pos);


fa
meandeltaAicc_30_zero = mean(AICc(:,1,1)) - mean(AICc(:,4,1)) ; 
meanevidenceratio_30_zero = exp((meandeltaAicc_30_zero)/2) ;
meandeltaAicc_5_zero = mean(AICc(:,1,2)) - mean(AICc(:,4,2)) ; 
meanevidenceratio_5_zero = exp((meandeltaAicc_5_zero)/2) ;
meandeltaAicc_30_large = mean(AICc(:,1,3)) - mean(AICc(:,4,3)) ; 
meanevidenceratio_30_large = exp((meandeltaAicc_30_large)/2); 
meandeltaAicc_5_large = mean(AICc(:,1,4)) - mean(AICc(:,4,4)) ; 
meanevidenceratio_5_large = exp((meandeltaAicc_5_large)/2) ;

result = mean([meanevidenceratio_30_zero, meanevidenceratio_5_zero, meanevidenceratio_30_large, meanevidenceratio_5_large])  


figure
for con = 1:4
    if con == 1
        C = C11; c = 0;
    elseif con == 2
        C = C22; c = 0.2;
    elseif con == 3
        C = C33; c = 0;
    elseif con == 4
        C = C44; c = 0.2;
    end    
    plot(1+c:nmodel+c, mean(R2adj(:,:,con)), '-o', 'Color', C, 'MarkerEdgeColor', C, 'MarkerFaceColor', C, 'linewidth',1,'MarkerSize',ms); hold on
%     errorbar(1+c:nmodel+c, mean(R2adj(:,:,con)), std(R2adj(:,:,con))/sqrt(N), 'o', 'color', C, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
    ylabel('R2adj');
    xticks(1+0.1:1:nmodel+0.1); xlim([0.7 nmodel+0.4]);
    set(gca,'XTickLabel',{'N(sigma)','L(b)','U(a,b)','Smax','S50','Smin'});
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.016 0]);
    title(tn(con));    
    legend({'30 dots, 0 penalty', '5 dots, 0 penalty', '30 dots, -500 penalty', '5 dots, -500 penalty'});        
end


