clear; close all;

load('dat_decision');

%%
for condition = 1:4    
    endpointorg = EndpointorgY(:,:,condition);
    samplesY = SamplesY(:,:,condition);
    samplenoise = SampleNoiseY(:,:,condition);
    ssy = []; sy = []; ssamples = [];
    
    [N T] = size(endpointorg); counts = [];
    for subi = 1:N
        for itrial = 1:T
            % each samples sorted in a order
            sy = samplesY(subi, itrial*31-30:itrial*31-1) ; %X1-X30
            endpoint = endpointorg(subi, :);

            if condition == 3 || condition == 4
                sy = sy(1:5);
            end
            isabovepenalty = sy > 180;
            counts(subi,itrial) = sum(isabovepenalty);
        end
    end
    meancounts(:,condition) = mean(counts,2);
    if condition == 1 || condition == 2
        pctgcounts(:,condition) = mean(counts,2) / 30 ;
    else
        pctgcounts(:,condition) = mean(counts,2) / 5 ;        
    end
end 
tmp = [1 3 2 4]; % [30,0], [5,0], [30,-500], [5,-500]
meancounts = meancounts(:,tmp);
ms = 5.5; cs = 10; ms2 = 10;
figure
C1 = [204 247 255]./255 ; C11 = [80 180 255]./255 ; 
C3 = [179 179 255]./255; C33 = [25 25 255]./255;
C2 = [255 221 200]./255; C22 = [255 140 0]./255;
C4 = [255 200 200]./255; C44 = [200 0 0]./255; 

ic = [1 2 3 4] + rand(17,1)*0.3-0.15; 
p = 0.3;
% plot(1:4, meancounts, 'k-'); hold on 
for j = 1:17 
    for i = 1:4
       plot(ic(j,:), meancounts(j,:), '-','Color',[0.5 0.5 0.5]); hold on 
    end
end
plot(ic(:,1), meancounts(:,1), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', C11, 'MarkerSize',ms); hold on
plot(ic(:,2), meancounts(:,2), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', C33, 'MarkerSize',ms); 
plot(ic(:,3), meancounts(:,3), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', C22, 'MarkerSize',ms); 
plot(ic(:,4), meancounts(:,4), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', C44, 'MarkerSize',ms); 

yticks(0:1:7); ylim([-0.25 6.5]); 
xlim([0.35 4.5]); xticks(1:1:5);
xticklabels({'30 dots, 0 penalty', '5 dots, 0 penalty', '30 dots, -500 penalty', '5 dots, -500 penalty'}); 
ylabel('The number of exception points in the penalty');
myfigAI2(13,0.025,400,430);


%% 
for condition = 1:4    
    endpointorg = EndpointorgY(:,:,condition);
    samplesY = SamplesY(:,:,condition);
    samplenoise = SampleNoiseY(:,:,condition);
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
        end
    end
    
    for subi = 1:N
        endpoint = endpointorg(subi, :)';
        sample = ssamples(:,:,subi);
        orgsample = sample - endpoint ;

        Smax(:, condition, subi) = orgsample(:,5) ;
        
        Endpoint(:, condition, subi) = endpoint ;        
    end
end 

% model fitting
for subi = 1:17

    Subdata   = Endpoint(:,:,subi) ;
    smax = Smax(:,:,subi) ;
    
    [estimates_smax, log_likelihood] = ML_fit_SmaxSt(Subdata, smax);
   
    Bias1(subi,:) = estimates_smax;
    Yhut1 = []; 
    
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
end

refboundary0 = Estimates_smax(:,1); % recovered reference boundary in 0 penalty condition
refboundary500 = Estimates_smax(:,2); % recovered reference boundary in -500 penalty condition

figure
plot(1:2,[refboundary0, refboundary500],'ko-')
mean([refboundary0, refboundary500])
std([refboundary0, refboundary500])
% unit: mm 




