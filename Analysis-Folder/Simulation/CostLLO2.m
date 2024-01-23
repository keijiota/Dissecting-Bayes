clear all; close all;

x = -10:0.05:5; %cm in steps of 0.5 mm
y = -5:0.05:5; 

penalty_rad = 2;  reward_rad = 2;

logodds = @(p) log(p./(1-p));
pm = [0.87 0.72];
Sigma = [1 0; 0 1]; %1sd~10mm=1cm
nrep = 100; 
simresult = [];

for is = 1:3
    if is == 1
       penalty_center = [-5 0]; reward_center = [-3 0];
       penalty = -500; reward = 100;
    elseif is == 2
       penalty_center = [-5 0]; reward_center = [-2 0];
       penalty = -500; reward = 100;        
    elseif is == 3
       penalty_center = [-5 0]; reward_center = [-1 0];
       penalty = -500; reward = 100;        
    elseif is == 4
       penalty_center = [-5 0]; reward_center = [-3 0];
       penalty = -100; reward = 100;
    elseif is == 5
       penalty_center = [-5 0]; reward_center = [-2 0];
       penalty = -100; reward = 100;        
    elseif is == 6
       penalty_center = [-5 0]; reward_center = [-1 0];
       penalty = -100; reward = 100;               
    end    

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

for imux = 1:length(x)
        imuy = round(length(y)/2) ;
        mux = x(imux); muy = y(imuy);
        
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
        eg_llo(imux, imuy) = prob_llo(1)*(penalty+reward) + prob_llo(2)*penalty + prob_llo(3)*reward;

end
[maxv a b] = mymax(eg);
[maxvv aa bb] = mymax(eg_llo); 
simresult = [simresult; maxv x(a) y(b) eg(aa,bb) maxvv x(aa) y(bb)];


th = 0:pi/100:2*pi;
sp = penalty_rad*cos(th)+ penalty_center(1);
cp = penalty_rad*sin(th)+ penalty_center(2);
sr = reward_rad*cos(th)+ reward_center(1);
cr = reward_rad*sin(th)+ reward_center(2);

subplot(1,3,is)
plot(sp,cp, 'r-'); hold on 
plot(sr,cr, 'b-'); 
plot(x(a),y(b),'k*');
plot(x(aa),y(bb),'gs');
axis('equal');
xlim([min(x) max(x)]); ylim([min(y) max(y)]);
box off;

end
pos(3) = 600; pos(4) = 600;
set(gcf, 'Position', pos);

simresult(:,4)./simresult(:,1)
simresult(:,5)./simresult(:,1)


