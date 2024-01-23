clear all; close all;

x = -20:0.05:25; %cm in steps of 0.5 mm
y = -15:0.05:15; 

penalty_rad = 8.97;  reward_rad = 8.97;

logodds = @(p) log(p./(1-p));
pm = [0.88 0.72];
s = mean([12, 23, 11, 20, 11].^(1/2)); 
Sigma = [s 0; 0 s]; %1sd~10mm=1cm
simresult = [];

penalty_center = [-2 0]; reward_center = [-2+penalty_rad*1.5 0];
penalty = -100; reward = 100;

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

imuy = round(length(y)/2);
muy = y(imuy);

nrep = 1000000; samples = ndRandn([0 0]', Sigma.^2, nrep)';
EG = [];
for ix = 10:0.01:15
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
    
    EG = [EG; ix, eg, eg_llo];    
end
plot(EG(:,1),EG(:,2)); hold on 
plot(EG(:,1),EG(:,3)); 
[a b] = mymax(EG(:,2)); optx  = EG(b,1);
[a b] = mymax(EG(:,3)); optlx = EG(b,1);


optx = 13.2; optlx = 13.3; 

nrep = 100; samples = ndRandn([optx 0]', Sigma.^2, nrep)';
mean(samples)

dis2penalty = sum((samples - penalty_center).^2, 2).^(1/2);
dis2reward  = sum((samples - reward_center).^2, 2).^(1/2);

prob_both = sum(dis2penalty <= penalty_rad & dis2reward <= reward_rad) / nrep ;
prob_penalty = sum(dis2penalty <= penalty_rad & dis2reward > reward_rad) / nrep ;
prob_reward = sum(dis2penalty > penalty_rad & dis2reward <= reward_rad) / nrep ;
prob_none = sum(dis2penalty > penalty_rad & dis2reward > reward_rad) / nrep ;
eg = prob_both*(penalty+reward) + prob_penalty*penalty + prob_reward*reward;

result(1,:) = [eg prob_both prob_penalty prob_reward prob_none];

xx = [prob_both prob_penalty prob_reward prob_none];
prob_llo = pm(1) * logodds(xx) + (1-pm(1)) * logodds(pm(2));
prob_llo = exp(prob_llo) ./ (1 + exp(prob_llo)) ;
eg_llo = prob_llo(1)*(penalty+reward) + prob_llo(2)*penalty + prob_llo(3)*reward;

result(2,:) = [eg_llo prob_llo]

th = 0:pi/100:2*pi;
sp = penalty_rad*cos(th)+ penalty_center(1);
cp = penalty_rad*sin(th)+ penalty_center(2);
sr = reward_rad*cos(th)+ reward_center(1);
cr = reward_rad*sin(th)+ reward_center(2);

figure
plot(sp,cp, 'r-'); hold on 
plot(sr,cr, 'b-'); 
plot(samples(:,1),samples(:,2),'k.', 'MarkerSize',11);
plot(optx, muy, 'rs', 'MarkerSize',6); 
plot(optlx, muy, 'bs', 'MarkerSize',6); 
axis('equal');
xlim([min(x) max(x)]); ylim([min(y) max(y)]);
box off;

pos(3) = 600; pos(4) = 600;
set(gcf, 'Position', pos);

% save s1 samples
