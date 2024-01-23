clear all; close all;

% [x,y] = meshgrid(-2:0.05:2);

x = -10:0.05:10;
y = -9:0.05:9;

penalty_rad = 2; penalty_center = [-5 0];
penalty = -500;
reward_rad = 2; reward_center = [-1 0];
reward = 100;

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

logodds = @(p) log(p./(1-p));
pm = [0.88 0.8];
Sigma = [1 0; 0 1];
nrep = 5000; yy = [-5:0.05:5];
for imux = 1:length(x)
    for imuy = 1:length(yy)
        mux = x(imux); muy = yy(imuy);
        
        samples = ndRandn([mux muy]', Sigma.^2, nrep)';
        dis2penalty = sum((samples - penalty_center).^2, 2).^(1/2);
        dis2reward  = sum((samples - reward_center).^2, 2).^(1/2);
        
        prob_both = sum(dis2penalty <= penalty_rad & dis2reward <= reward_rad) / nrep ;
        prob_penalty = sum(dis2penalty <= penalty_rad & dis2reward > reward_rad) / nrep ;
        prob_reward = sum(dis2penalty > penalty_rad & dis2reward <= reward_rad) / nrep ;
        prob_none = sum(dis2penalty > penalty_rad & dis2reward > reward_rad) / nrep ;
        eg(imux, imuy) = prob_both*(penalty+reward) + prob_penalty*penalty + prob_reward*reward;
        
        xx = [prob_both prob_penalty prob_reward prob_none];
        prob_llo = pm(1) * logodds(xx) + (1-pm(1)) * logodds(pm(2));
        prob_llo = exp(prob_llo) ./ (1 + exp(prob_llo)) ;
%         prob_llo = [prob_llo(1) prob_llo(2) prob_reward prob_none];
        eg_llo(imux, imuy) = prob_llo(1)*(penalty+reward) + prob_llo(2)*penalty + prob_llo(3)*reward;

    end
end
[maxv a b]= mymax(eg);
[maxvv aa bb]= mymax(eg_llo);
[maxv x(a) yy(b)]
[eg(aa,bb) maxvv x(aa) yy(bb)]

% scatter(samples(:,1),samples(:,2));
% xlim([min(x) max(x)]); ylim([min(y) max(y)]);

figure
subplot(1,3,1)
mesh(x,y,Gainfunc');  view(0,90); colorbar; ylim([min(y) max(y)]);
subplot(1,3,2)
mesh(x,yy,eg');  view(0,90); colorbar; ylim([min(y) max(y)]);
hold on
plot3(x(a),yy(b),maxv,'rs');
subplot(1,3,3)
mesh(x,yy,eg_llo');  view(0,90); colorbar; ylim([min(y) max(y)]);
hold on
plot3(x(aa),yy(bb),maxvv,'rs');

