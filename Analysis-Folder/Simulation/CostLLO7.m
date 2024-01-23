clear all; close all;

x = -20:0.05:25; %cm in steps of 0.5 mm
y = -15:0.05:15; 

penalty_rad = 8.97;  reward_rad = 8.97;
s = mean([11.94, 23.36, 11.07, 19.65, 11.43].^(1/2));
Sigma = [s 0; 0 s]; %1sd~10mm=1cm

penalty_center = [-2 0]; reward_center = [-2+penalty_rad*1.5 0]; % 1, 1.5, or 2R
penalty = -100; % 0, -100, -500 
reward = 100; % 0, 100, 500

Gainfunc = []; eg = []; eg_llo = [];
for ix = 1:length(x)
    for iy = 1:length(y)
        xy = [x(ix) y(iy)];
        dis2penalty = sum((xy - penalty_center).^2)^(1/2);
        dis2reward = sum((xy - reward_center).^2)^(1/2);
        
        if dis2penalty <= penalty_rad && dis2reward <= reward_rad
            Gainfunc(ix,iy) = penalty + reward;
        elseif dis2penalty < penalty_rad
            Gainfunc(ix,iy) = penalty ;
        elseif dis2reward < reward_rad
            Gainfunc(ix,iy) = reward;
        else
            Gainfunc(ix,iy) = 0;
        end
    end
end

imuy = round(length(y)/2); muy = y(imuy);

nrep = 10000000; samples = ndRandn([0 0]', Sigma.^2, nrep)';
logodds = @(p) log(p./(1-p));
pm = [0.88 0.76];

% EG from pdf (~= many emprical samples)

% EG = [];
% for ix = -5:0.05:20
%     samples_s = [samples(:,1)+ix samples(:,2)]; 
%     dis2penalty = sum((samples_s - penalty_center).^2, 2).^(1/2); dis2reward = sum((samples_s - reward_center).^2, 2).^(1/2);
%     
%     prob_both = sum(dis2penalty <= penalty_rad & dis2reward <= reward_rad) / nrep ;
%     prob_penalty = sum(dis2penalty <= penalty_rad & dis2reward > reward_rad) / nrep ;
%     prob_reward = sum(dis2penalty > penalty_rad & dis2reward <= reward_rad) / nrep ;
%     prob_none = sum(dis2penalty > penalty_rad & dis2reward > reward_rad) / nrep ;
%     eg = prob_both*(penalty+reward) + prob_penalty*penalty + prob_reward*reward;
%     
%     xx = [prob_both prob_penalty prob_reward prob_none];
%     prob_llo = pm(1) * logodds(xx) + (1-pm(1)) * logodds(pm(2));
%     prob_llo = exp(prob_llo) ./ (1 + exp(prob_llo)) ;
%     eg_llo = prob_llo(1)*(penalty+reward) + prob_llo(2)*penalty + prob_llo(3)*reward;
%     
%     EG = [EG; ix, eg, eg_llo, prob_both, prob_penalty, prob_reward, prob_none, prob_llo];    
% end
% plot(EG(:,1),EG(:,2)); hold on 
% plot(EG(:,1),EG(:,3)); 
% [opteg tmp1] = mymax(EG(:,2)); optx = EG(tmp1,1); optx = median(optx); 
% tmp = find(abs(EG(:,1)-round(optx,1))<0.01); prob = EG(tmp,4:end);
% 
% [optegllo tmp2] = mymax(EG(:,3)); optlx = EG(tmp2,1); optlx = median(optlx);
% tmp = find(abs(EG(:,1)-round(optlx,1))<0.01); probllo = EG(tmp,4:end);
% result = [opteg optx prob; optegllo optlx probllo]


optx = 13.15; optlx = 13.30; % by simulation 
EG = [];
for ix = [optx optlx]
    samples_s = [samples(:,1)+ix samples(:,2)]; 
    dis2penalty = sum((samples_s - penalty_center).^2, 2).^(1/2); dis2reward = sum((samples_s - reward_center).^2, 2).^(1/2);
    
    prob_both = sum(dis2penalty <= penalty_rad & dis2reward <= reward_rad) / nrep ;
    prob_penalty = sum(dis2penalty <= penalty_rad & dis2reward > reward_rad) / nrep ;
    prob_reward = sum(dis2penalty > penalty_rad & dis2reward <= reward_rad) / nrep ;
    prob_none = sum(dis2penalty > penalty_rad & dis2reward > reward_rad) / nrep ;
    eg = prob_both*(penalty+reward) + prob_penalty*penalty + prob_reward*reward;
    
    xx = [prob_both prob_penalty prob_reward prob_none];
    prob_llo = pm(1) * logodds(xx) + (1-pm(1)) * logodds(pm(2));
    prob_llo = exp(prob_llo) ./ (1 + exp(prob_llo)) ;
    eg_llo = prob_llo(1)*(penalty+reward) + prob_llo(2)*penalty + prob_llo(3)*reward;
    
    EG = [EG; ix, eg, eg_llo, prob_both, prob_penalty, prob_reward, prob_none, prob_llo];    
end

xx = [0.035 0.002 0.879 0.084]; % sum to be 1
prob_llo = pm(1) * logodds(xx) + (1-pm(1)) * logodds(pm(2));
prob_llo = exp(prob_llo) ./ (1 + exp(prob_llo)) ;

[EG(1,2) xx; EG(1,3) prob_llo]


% endpoint
nrep = 100; samples = ndRandn([0 0]', Sigma.^2, nrep)';
load s4; 

sampless = [samples(:,1)+optx, samples(:,2)];
    
th = 0:pi/100:2*pi;
sp = penalty_rad*cos(th)+ penalty_center(1);
cp = penalty_rad*sin(th)+ penalty_center(2);
sr = reward_rad*cos(th)+ reward_center(1);
cr = reward_rad*sin(th)+ reward_center(2);

figure
plot(sp,cp, 'r-'); hold on 
plot(sr,cr, 'b-'); 
plot(sampless(:,1),sampless(:,2),'k.', 'MarkerSize',11);
plot(optx, muy, 'rs', 'MarkerSize',6); 
plot(optlx, muy, 'bs', 'MarkerSize',6); 
axis('equal');
xlim([min(x) max(x)]); ylim([min(y) max(y)]);
box off;
pos(3) = 600; pos(4) = 600;
set(gcf, 'Position', pos);
% save s4 samples

% 
figure; hold on
x = 0:0.01:1;
y = pm(1) * logodds(x) + (1-pm(1)) * logodds(pm(2));
y = exp(y) ./ (1 + exp(y)) ;
plot(xx,prob_llo, 'ko', 'MarkerFaceColor','w', 'MarkerSize',7) ; hold on
plot(0:0.1:1.2,0:0.1:1.2,'k--');
plot(x, y, 'k-');
axis('square');
xlim([0 1.05]); ylim([0 1.05]);
xticks(0:0.2:1); yticks(0:0.2:1);
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);
box off; 
pos(3) = 200; pos(4) = 200;
set(gcf, 'Position', pos);











