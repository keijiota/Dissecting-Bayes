clear; close all ;

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

%% making a matrix we need
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
            
            mu_hut = mean(sy);
            Sigma_hut(subi,itrial) = std(sy);
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
        % given sample variance and N, the posterior variance of estimated population variance can be defined 
        % given that estimated population variance, the estimated pdfs can be defiend

    end
end 
load optimalaimpoint

tn = {'30 dots & 0','30 dots & -500','5 dots & 0','5 dots & -500'};

%% model fitting
for subi = 1:N
    Optaim_sg  = reshape(OptAim_Sufficient_G(subi,:,:),50,4) ; 
    Subdata   = Endpoint(:,:,subi) ;
    smax = Smax(:,:,subi) ;
    s50 = S50(:,:,subi) ;
    smin = Smin(:,:,subi) ;
    
    [estimates_smax, log_likelihood] = ML_fit_MaxPoint(Subdata, smax);
   
    Bias1(subi,:) = estimates_smax;
    Yhut1 = []; Yhut2 = []; Yhut3 = [];    
    
%     figure
    for con = 1:4
        if con == 1 || con == 3
            bias1 = estimates_smax(1);
        else
            bias1 = estimates_smax(2);
        end
        Yhut1(:,con) = (180 - bias1) - smax(:,con) ;
    end
    
    YHut1(:,:,subi) = Yhut1;
    Estimates_smax(subi,:) = estimates_smax;
    
    for con = 1:4
        subdata = Subdata(:,con);
        optaim_sg = Optaim_sg(:,con); 
        yhut1 = Yhut1(:,con);
        [T M] = size(subdata);

        % argmax G(x)*P(x)*P(PopSigma|SampleSigma)
        diff = sum((subdata - optaim_sg).^2);
        log_likelihood = (T*M)/2 * (log(2*pi) + log(diff) + 1) ;
        K = 0; c = 1;
        AICc(subi, c, con) = -2*(-log_likelihood) + 2*K + 2*(K*(K+1))/(T*M-K-1) ;
        BIC(subi, c, con) = -2*(-log_likelihood) + log(T*M)*K;
        R2(subi, c, con) = 1 - (sum((optaim_sg - subdata).^2) / sum((mean(subdata) - subdata).^2)) ;
        R2adj(subi, c, con) = 1 - ((1-R2(subi, c, con))*(T*M - 1) / (T*M - K - 1)) ;                  
                
        % y = (PB - bias_gf) - Max point strategy
        diff = sum((subdata - yhut1).^2);
        log_likelihood = (T*M)/2 * (log(2*pi) + log(diff) + 1) ;
        K = 1; c = 2;
        AICc(subi, c, con) = -2*(-log_likelihood) + 2*K + 2*(K*(K+1))/(T*M-K-1) ;
        BIC(subi, c, con) = -2*(-log_likelihood) + log(T*M)*K;               
        R2(subi, c, con) = 1 - (sum((yhut1 - subdata).^2) / sum((mean(subdata) - subdata).^2)) ;
        R2adj(subi, c, con) = 1 - ((1-R2(subi, c, con))*(T*M - 1) / (T*M - K - 1)) ;
        
    end
end

for subi = 1:N
    avSubdata(subi,:) = mean(Endpoint(:,:,subi));    
    avOptaim_sg(subi,:) = squeeze(mean(OptAim_Sufficient_G(subi,:,:)));        
    avYHut1(subi,:) = mean(YHut1(:,:,subi));
end

tmp = [1 3 2 4]; % [30,0], [5,0], [30,-500], [5,-500]
avSubdata = avSubdata(:,tmp);
avOptaim_sg = avOptaim_sg(:,tmp);
avYHut1 = avYHut1(:,tmp);
AICc = AICc(:,:,tmp) ; 
BIC = BIC(:,:,tmp) ; 
R2adj = R2adj(:,:,tmp) ; 

%% 
C1 = [204 247 255]./255 ; C11 = [100 180 255]./255 ; 
C2 = [179 179 255]./255; C22 = [70 70 255]./255;
C3 = [255 221 150]./255; C33 = [255 160 20]./255;
C4 = [255 200 200]./255; C44 = [220 60 60]./255; 
C5 = [100 100 100]./255;
cs = 12; ms = 9.5; ms2 = 11;

figure
%--------Supplementary Figure 7A---------------
ic = [1 2.2 3.4 4.6]; p = 0.3;
for con = 1:4
    if con == 1
        C = C11; CC = C11; 
    elseif con == 2
        C = C22; CC = C22; 
    elseif con == 3
        C = C33; CC = C33; 
    elseif con == 4
        C = C44; CC = C44; 
    end
    subplot(1,58,1:10.5); ip = 14;
    plot(ic(con), avSubdata(:,con), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', C, 'MarkerSize',ms,'linewidth',0.1); hold on
    plot(ic(con)+p, mean(avSubdata(:,con)), 's', 'MarkerEdgeColor', CC, 'MarkerFaceColor', CC, 'MarkerSize',ms2); hold on
    errorbar(ic(con)+p, mean(avSubdata(:,con)), 2*std(avSubdata(:,con))/sqrt(17), 'o', 'color', CC, 'linewidth', 1.5, 'MarkerSize', 0.1, 'CapSize',cs);
    box off;
end
plot(ic(1:2)+p, [mean(avSubdata(:,1)), mean(avSubdata(:,2))], 'k--', 'linewidth', 1);
plot(ic(3:4)+p, [mean(avSubdata(:,3)), mean(avSubdata(:,4))], 'k--', 'linewidth', 1);
plot(ic(1:2)+p, [mean(avOptaim_sg(:,1)), mean(avOptaim_sg(:,2))], 'dk--', 'linewidth', 1,'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'k', 'MarkerSize',ms+1);
plot(ic(3:4)+p, [mean(avOptaim_sg(:,3)), mean(avOptaim_sg(:,4))], 'dk--', 'linewidth', 1,'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'k', 'MarkerSize',ms+1);

ylim([110 180]); xlim([ic(1)-0.75 ic(4)+0.75]); xticks(ic+p/2); yticks(100:20:180);
xticklabels({'30, 0', '5, 0', '30, -500', '5, -500'});
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);


%--------Figure 9 C&D---------------
ms = 9;
for iplot = 1:2
    if iplot == 1
        dat = avOptaim_sg; yn = 'BDT model';
    elseif iplot == 2
        dat = avYHut1; yn = 'Max-point model';
    end
    
    for con = 1:4
        if con == 1
            C = C11; 
        elseif con == 2
            C = C22;  
        elseif con == 3
            C = C33;  
        elseif con == 4
            C = C44;  
        end
        subplot(1,58,ip:ip+12);
        plot(dat(:,con), avSubdata(:,con),  'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', C, 'MarkerSize',ms,'linewidth',0.01); hold on
    end
    
    ip = ip+12+4;
    plot(-110:180,-110:180,'k--');
    ylim([110 180]); xlim([110 180]); xticks(120:20:180); yticks(120:20:180);
    axis square
    ylabel('Actual set point');
    xlabel(yn, 'FontName', 'Arial', 'FontSize', 10);
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);
    box off;    
end
legend({'30 dots, 0 penalty','5 dots, 0 penalty','30 dots, -500 penalty','5 dots, -500 penalty'})

%--------Figure 9B---------------
figure; 
tmp = [1 2]; 
nmodel = length(tmp);
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
    plot(1+c:nmodel+c, mean(AICc(:,tmp,con)), '-s', 'Color', C, 'MarkerEdgeColor', C, 'MarkerFaceColor', C, 'linewidth',1,'MarkerSize',ms2); hold on
end
ylabel('AICc');
xticks(1+0.15:1:nmodel+0.15);  xlim([0.6 2.6]); ylim([540 620]); yticks(540:20:640);
set(gca,'XTickLabel',{'BDT','Max-point'});
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.016 0]);
% pos(3) = 300; pos(4) = 400; set(gcf, 'Position', pos);

avdiffAICc = mean(mean(AICc(:,tmp,:)),3)
avdiffBIC = mean(mean(BIC(:,tmp,:)),3)

% Sufficient G vs Smax
meandeltaAicc_30_zero = mean(AICc(:,1,1)) - mean(AICc(:,2,1)) ; 
meanevidenceratio_30_zero = exp((meandeltaAicc_30_zero)/2) ;
meandeltaAicc_5_zero = mean(AICc(:,1,2)) - mean(AICc(:,2,2)) ; 
meanevidenceratio_5_zero = exp((meandeltaAicc_5_zero)/2) ;
meandeltaAicc_30_large = mean(AICc(:,1,3)) - mean(AICc(:,2,3)) ; 
meanevidenceratio_30_large = exp((meandeltaAicc_30_large)/2); 
meandeltaAicc_5_large = mean(AICc(:,1,4)) - mean(AICc(:,2,4)) ; 
meanevidenceratio_5_large = exp((meandeltaAicc_5_large)/2) ;

result = mean([meanevidenceratio_30_zero, meanevidenceratio_5_zero, meanevidenceratio_30_large, meanevidenceratio_5_large])  




