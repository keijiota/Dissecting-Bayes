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

        % Gaussian
        %          endpoint_dist_hut = exp(-(y - Y(iy) + 0).^2./(2*v_hut')) ./ sqrt(2*pi*v_hut');
        
        % Laplacian: variance of Laplacian distribution ->  sigma = 2*b^2  =>  b = sqrt(sigma)/sqrt(2)
        %         b_hut = sqrt(v_hut)/sqrt(2);
        %         endpoint_dist_hut_L = 1./(2*b_hut') .* exp(-(abs(y-Y(iy))./b_hut')) ;
        
        % Uniform: variance of uniform distribution ->  sigma = (b-a)^2 / 12  =>  b-a = 2*sqrt(3)*sqrt(sigma)
        %         b_a = 2*sqrt(3)*sqrt(v_hut);
        %         bu = Y(iy) + (b_a'/2); au = Y(iy) - (b_a'/2);
        %         endpoint_dist_hut_U = zeros(length(v_hut),length(y));
        %         tmp = y>au & y<bu;
        %         for iv = 1:length(v_hut)
        %             endpoint_dist_hut_U(iv, tmp(iv,:)) = 1/(bu(iv) - au(iv));
        %         end
    end
end 

load suffdata

tn = {'30 dots & 0','30 dots & -500','5 dots & 0','5 dots & -500'};

%% model fitting
for subi = 1:N
    Optaim_sg  = reshape(OptAim_Sufficient_G(subi,:,:),50,4) ; 
    Optaim_sl  = reshape(OptAim_Sufficient_L(subi,:,:),50,4) ; 
    Optaim_su  = reshape(OptAim_Sufficient_U(subi,:,:),50,4) ; 
    Subdata   = Endpoint(:,:,subi) ;
    smax = Smax(:,:,subi) ;
    s50 = S50(:,:,subi) ;
    smin = Smin(:,:,subi) ;
    
    [estimates_smax, log_likelihood] = ML_fit_SmaxSt(Subdata, smax);
    [estimates_s50,  log_likelihood] = ML_fit_S50St(Subdata, s50);
    [estimates_smin, log_likelihood] = ML_fit_SminSt(Subdata, smin);
   
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
    end
    
    YHut1(:,:,subi) = Yhut1;
    YHut2(:,:,subi) = Yhut2;
    YHut3(:,:,subi) = Yhut3;    
    Estimates_smax(subi,:) = estimates_smax;
    Estimates_s50(subi,:) = estimates_s50;    
    Estimates_smin(subi,:) = estimates_smin;
    
    for con = 1:4
        subdata = Subdata(:,con);
        optaim_sg = Optaim_sg(:,con);        
        optaim_sl = Optaim_sl(:,con);
        optaim_su = Optaim_su(:,con);
        yhut1 = Yhut1(:,con);
        yhut2 = Yhut2(:,con);
        yhut3 = Yhut3(:,con);
        [T M] = size(subdata);

        % argmax G(x)*P(x)*P(PopSigma|SampleSigma)
        diff = sum((subdata - optaim_sg).^2);
        log_likelihood = (T*M)/2 * (log(2*pi) + log(diff) + 1) ;
        K = 0; c = 1;
        AICc(subi, c, con) = -2*(-log_likelihood) + 2*K + 2*(K*(K+1))/(T*M-K-1) ;
        BIC(subi, c, con) = -2*(-log_likelihood) + log(T*M)*K;
        R2(subi, c, con) = 1 - (sum((optaim_sg - subdata).^2) / sum((mean(subdata) - subdata).^2)) ;
        R2adj(subi, c, con) = 1 - ((1-R2(subi, c, con))*(T*M - 1) / (T*M - K - 1)) ;                  
        
        %  argmax G(x)*P(x)*P(PopSigma|SampleSigma), P(x)~Laplacian
        diff = sum((subdata - optaim_sl).^2);
        log_likelihood = (T*M)/2 * (log(2*pi) + log(diff) + 1) ;
        K = 0; c = 2;
        AICc(subi, c, con) = -2*(-log_likelihood) + 2*K + 2*(K*(K+1))/(T*M-K-1) ;
        BIC(subi, c, con) = -2*(-log_likelihood) + log(T*M)*K;
        R2(subi, c, con) = 1 - (sum((optaim_sl - subdata).^2) / sum((mean(subdata) - subdata).^2)) ;
        R2adj(subi, c, con) = 1 - ((1-R2(subi, c, con))*(T*M - 1) / (T*M - K - 1)) ;
        
        %  argmax G(x)*P(x)*P(PopSigma|SampleSigma), P(x)~Uniform
        diff = sum((subdata - optaim_su).^2);
        log_likelihood = (T*M)/2 * (log(2*pi) + log(diff) + 1) ;
        K = 0; c = 3;
        AICc(subi, c, con) = -2*(-log_likelihood) + 2*K + 2*(K*(K+1))/(T*M-K-1) ;
        BIC(subi, c, con) = -2*(-log_likelihood) + log(T*M)*K;                       
        R2(subi, c, con) = 1 - (sum((optaim_su - subdata).^2) / sum((mean(subdata) - subdata).^2)) ;
        R2adj(subi, c, con) = 1 - ((1-R2(subi, c, con))*(T*M - 1) / (T*M - K - 1)) ;
        
        % y = (PB - bias_gf) - Smax
        diff = sum((subdata - yhut1).^2);
        log_likelihood = (T*M)/2 * (log(2*pi) + log(diff) + 1) ;
        K = 1; c = 4;
        AICc(subi, c, con) = -2*(-log_likelihood) + 2*K + 2*(K*(K+1))/(T*M-K-1) ;
        BIC(subi, c, con) = -2*(-log_likelihood) + log(T*M)*K;               
        R2(subi, c, con) = 1 - (sum((yhut1 - subdata).^2) / sum((mean(subdata) - subdata).^2)) ;
        R2adj(subi, c, con) = 1 - ((1-R2(subi, c, con))*(T*M - 1) / (T*M - K - 1)) ;
        
        % y = (PB - bias_gf) - S50
        diff = sum((subdata - yhut2).^2);
        log_likelihood = (T*M)/2 * (log(2*pi) + log(diff) + 1) ;
        K = 1; c = 5;
        AICc(subi, c, con) = -2*(-log_likelihood) + 2*K + 2*(K*(K+1))/(T*M-K-1) ;
        BIC(subi, c, con) = -2*(-log_likelihood) + log(T*M)*K;                       
        R2(subi, c, con) = 1 - (sum((yhut2 - subdata).^2) / sum((mean(subdata) - subdata).^2)) ;
        R2adj(subi, c, con) = 1 - ((1-R2(subi, c, con))*(T*M - 1) / (T*M - K - 1)) ;
        
        % y = (PB - bias_gf) - Smin
        diff = sum((subdata - yhut3).^2);
        log_likelihood = (T*M)/2 * (log(2*pi) + log(diff) + 1) ;
        K = 1; c = 6;
        AICc(subi, c, con) = -2*(-log_likelihood) + 2*K + 2*(K*(K+1))/(T*M-K-1) ;
        BIC(subi, c, con) = -2*(-log_likelihood) + log(T*M)*K;                       
        R2(subi, c, con) = 1 - (sum((yhut3 - subdata).^2) / sum((mean(subdata) - subdata).^2)) ;
        R2adj(subi, c, con) = 1 - ((1-R2(subi, c, con))*(T*M - 1) / (T*M - K - 1)) ;
    end
end

for subi = 1:N
    avSubdata(subi,:) = mean(Endpoint(:,:,subi));    
    avOptaim_sg(subi,:) = squeeze(mean(OptAim_Sufficient_G(subi,:,:)));        
    avOptaim_sl(subi,:) = squeeze(mean(OptAim_Sufficient_L(subi,:,:)));        
    avOptaim_su(subi,:) = squeeze(mean(OptAim_Sufficient_U(subi,:,:)));        
    avYHut1(subi,:) = mean(YHut1(:,:,subi));
    avYHut2(subi,:) = mean(YHut2(:,:,subi));
    avYHut3(subi,:) = mean(YHut3(:,:,subi));        
end

tmp = [1 3 2 4]; % [30,0], [5,0], [30,-500], [5,-500]
avSubdata = avSubdata(:,tmp);
avOptaim_sg = avOptaim_sg(:,tmp);
avOptaim_sl = avOptaim_sl(:,tmp);
avOptaim_su = avOptaim_su(:,tmp);
avYHut1 = avYHut1(:,tmp);
avYHut2 = avYHut2(:,tmp);
avYHut3 = avYHut3(:,tmp);
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

ms = 9;
for iplot = 1:3
    if iplot == 1
        dat = avOptaim_sg; yn = 'Sufficient model';
    elseif iplot == 2
        dat = avYHut1; yn = 'Smax model';
    elseif iplot == 3
        dat = avYHut3; yn =  'Smin model';
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


figure
tmp = [1:4]; % leaving S50%, Smin
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
xticks(1+0.15:1:nmodel+0.15); xlim([0.25 nmodel+0.75]); ylim([540 630]); yticks(540:20:630);
set(gca,'XTickLabel',{'N(Popsigma)','L(b)','U(a,b)','Smax'});
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.016 0]);
pos(3) = 600; pos(4) = 400;
set(gcf, 'Position', pos);

figure; subplot(121);
tmp = [1 4]; 
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
set(gca,'XTickLabel',{'N(Popsigma)','Smax'});
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.016 0]);
% pos(3) = 300; pos(4) = 400; set(gcf, 'Position', pos);

subplot(122);
tmp = [1 4]; 
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
    plot(1+c:nmodel+c, mean(BIC(:,tmp,con)), '-s', 'Color', C, 'MarkerEdgeColor', C, 'MarkerFaceColor', C, 'linewidth',1,'MarkerSize',ms2); hold on
end
ylabel('BIC');
xticks(1+0.15:1:nmodel+0.15); xlim([0.6 2.6]); ylim([540 620]); yticks(540:20:640);
set(gca,'XTickLabel',{'N(Popsigma)','Smax'});
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.016 0]);
pos(3) = 500; pos(4) = 400; set(gcf, 'Position', pos);

avdiffAICc = mean(mean(AICc(:,tmp,:)),3)
avdiffBIC = mean(mean(BIC(:,tmp,:)),3)



% Sufficient G vs Smax
meandeltaAicc_30_zero = mean(AICc(:,1,1)) - mean(AICc(:,4,1)) ; 
meanevidenceratio_30_zero = exp((meandeltaAicc_30_zero)/2) ;
meandeltaAicc_5_zero = mean(AICc(:,1,2)) - mean(AICc(:,4,2)) ; 
meanevidenceratio_5_zero = exp((meandeltaAicc_5_zero)/2) ;
meandeltaAicc_30_large = mean(AICc(:,1,3)) - mean(AICc(:,4,3)) ; 
meanevidenceratio_30_large = exp((meandeltaAicc_30_large)/2); 
meandeltaAicc_5_large = mean(AICc(:,1,4)) - mean(AICc(:,4,4)) ; 
meanevidenceratio_5_large = exp((meandeltaAicc_5_large)/2) ;

result = mean([meanevidenceratio_30_zero, meanevidenceratio_5_zero, meanevidenceratio_30_large, meanevidenceratio_5_large])  

% Smax G vs Smin
meandeltaAicc_30_zero = mean(AICc(:,6,1)) - mean(AICc(:,4,1)) ; 
meanevidenceratio_30_zero = exp((meandeltaAicc_30_zero)/2) ;
meandeltaAicc_5_zero = mean(AICc(:,6,2)) - mean(AICc(:,4,2)) ; 
meanevidenceratio_5_zero = exp((meandeltaAicc_5_zero)/2) ;
meandeltaAicc_30_large = mean(AICc(:,6,3)) - mean(AICc(:,4,3)) ; 
meanevidenceratio_30_large = exp((meandeltaAicc_30_large)/2); 
meandeltaAicc_5_large = mean(AICc(:,6,4)) - mean(AICc(:,4,4)) ; 
meanevidenceratio_5_large = exp((meandeltaAicc_5_large)/2) ;

result = mean([meanevidenceratio_30_zero, meanevidenceratio_5_zero, meanevidenceratio_30_large, meanevidenceratio_5_large])  


Result = [];
for i = 1:4
    [tvalue, df, p, h, d, meandiff, CI_l, CI_u] = myttest(avSubdata(:,i),avOptaim_sg(:,i));
    Result = [Result; tvalue, df, p, h, d, meandiff, CI_l, CI_u];    
end

Result2 = [];
A = avSubdata(:,2) - avSubdata(:,1);
B = avOptaim_sg(:,2) - avOptaim_sg(:,1);
[tvalue, df, p, h, d, meandiff, CI_l, CI_u] = myttest(A,B);
Result2 = [Result2; tvalue, df, p, h, d, meandiff, CI_l, CI_u];    

A = avSubdata(:,4) - avSubdata(:,3);
B = avOptaim_sg(:,4) - avOptaim_sg(:,3);
[tvalue, df, p, h, d, meandiff, CI_l, CI_u] = myttest(A,B);
Result2 = [Result2; tvalue, df, p, h, d, meandiff, CI_l, CI_u];    













