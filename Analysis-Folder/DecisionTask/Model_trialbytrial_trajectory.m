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
        
        for itrial = 1:T
            sigma_hut = sigmahut(itrial)*10; b_hut = bhut(itrial)*10;
            l_hut = lhut(itrial)*10; m_hut = mhut(itrial)*10;            
            EG_g = []; EG_l = []; EG_u = [];
            for mu = 1:1:2000
                gaussian = exp(-(x - mu).^2/(2*sigma_hut^2)) / sqrt(2*pi*sigma_hut^2) ;
                gaussian = gaussian / sum(gaussian);
                EG_g = [EG_g; sum(gaussian.*gainfunc), mu];
                
                laplacian = 1/(2*b_hut) * exp(-(abs(x-mu)/b_hut)) ;
                laplacian = laplacian / sum(laplacian);
                EG_l = [EG_l; sum(laplacian.*gainfunc), mu];
                
                Mu = (l_hut + m_hut)/2 + mu;
                uniform = zeros(1,length(x)) ;
                uniform(round(l_hut+mu)+abs(min(x)):round(m_hut+mu)+abs(min(x))) = 1/(l_hut - m_hut);
                uniform = uniform/sum(uniform) ;
                EG_u = [EG_u; sum(uniform.*gainfunc), Mu];
            end
            
            xx = 1:1:2000;
            [OptEG_g d] = max(EG_g(:,1));
            OptAim_g(itrial, condition, subi) = EG_g(d,2)/10 ;
            
            [OptEG_l d] = max(EG_l(:,1));
            OptAim_l(itrial, condition, subi) = EG_l(d,2)/10 ;
            
            [OptEG_u d] = max(EG_u(:,1));
            OptShift = xx(d);
            OptAim_u(itrial, condition, subi) = EG_u(d,2)/10 ;
        end        

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
    
    figure
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
        
        subplot(2,2,con)
        plot(Subdata(:,con), 'k-','linewidth',2); hold on
        plot(Optaim_g(:,con), '-','linewidth',0.5,'Color',[1 0 0]); 
        plot(Optaim_l(:,con), '-','linewidth',0.5, 'Color',[0 1 0]); 
        plot(Optaim_u(:,con), '-','linewidth',0.5, 'Color',[0 0 1]);         
        plot(Yhut1(:,con), '--','linewidth',1.5,'Color',[1 0 0]);
        plot(Yhut2(:,con), '--','linewidth',1.5,'Color',[0 1 0]);        
        plot(Yhut3(:,con), '--','linewidth',0.5, 'Color',[0 0 1]);         
        title(tn(con));
        if con == 4
           legend({'Data','N(sigma)','L(b)','U(a,b)','Smax','S50','Smin'}); 
        end
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

for con = 1:4
    trjavSubdata(con,:) = mean(squeeze(Endpoint(:,con,:))');     
    trjavOptaim_g(con,:) = mean(squeeze(OptAim_g(:,con,:))'); 
    trjavOptaim_l(con,:) = mean(squeeze(OptAim_l(:,con,:))'); 
    trjavOptaim_u(con,:) = mean(squeeze(OptAim_u(:,con,:))');     
    trjavYHut1(con,:) = mean(squeeze(YHut1(:,con,:))'); 
    trjavYHut2(con,:) = mean(squeeze(YHut2(:,con,:))'); 
    trjavYHut3(con,:) = mean(squeeze(YHut3(:,con,:))');     
end

Endpoint = [Endpoint(:,1,:), Endpoint(:,3,:), Endpoint(:,2,:), Endpoint(:,4,:)];
OptAim_g = [OptAim_g(:,1,:), OptAim_g(:,3,:), OptAim_g(:,2,:), OptAim_g(:,4,:)];
YHut1 = [YHut1(:,1,:), YHut1(:,3,:), YHut1(:,2,:), YHut1(:,4,:)];


lw = 1; alphaa = 0.4; 
figure
for con = 1:4
    subplot(2,2,con)
    seshade(squeeze(Endpoint(:,con,:))',alphaa,'k', 'k',1:50); hold on 
    seshade(squeeze(OptAim_g(:,con,:))',alphaa,'b', 'b',1:50); 
    seshade(squeeze(YHut1(:,con,:))',alphaa,'r', 'r',1:50); hold on 

%      plot(trjavSubdata(con,:), 'k-','linewidth',lw); hold on
%     plot(trjavOptaim_g(con,:), '-','linewidth',lw, 'Color',[0 0 1]);
%     plot(trjavOptaim_l(con,:), '--','linewidth',lw, 'Color',[0 0 1]);
%     plot(trjavOptaim_u(con,:), ':','linewidth',lw, 'Color',[0 0 1]);
%     plot(trjavYHut1(con,:), '-','linewidth',lw, 'Color',[1 0 0]);
%     plot(trjavYHut2(con,:), '--','linewidth',lw, 'Color',[1 0 0]);
%     plot(trjavYHut3(con,:), ':','linewidth',lw, 'Color',[1 0 0]);       

    ylabel('Set point'); xlabel('Trials'); 
    xticks(0:25:50); xlim([0 51]);
    ylim([120 180]); yticks(120:10:180);
    set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.016 0]);    
end




