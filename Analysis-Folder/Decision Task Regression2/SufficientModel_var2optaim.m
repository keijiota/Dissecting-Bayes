clear ; close all ;

y = -500:0.1:750;
% v_hut = [0:0.001:0.1, 0.2:0.1:500, 501:1:20000];
v_hut = [0:0.1:5000];
m_hut = -10:0.1:10; 
Y = 120:0.01:180;

TL = 180 ;
for i = 1:length(y)
    if y(i) < 0
        gainFunc(i) = 0;
        gainFunc2(i) = 0;
    elseif y(i) >= 0 && y(i) <= TL
        gainFunc(i) = y(i)/(TL)*100;
        gainFunc2(i) = y(i)/(TL)*100;
    elseif y(i) > TL
        gainFunc(i) = 0;
        gainFunc2(i) = -500;
    end
end

% prior probability
population_mean = [0 0];
population_var  = [10 20].^2;

inrange =  m_hut >= population_mean(1) & m_hut <= population_mean(end) ;
prior_mean = inrange * 1/sum(inrange);

inrange =  v_hut >= population_var(1) & v_hut <= population_var(end) ;
prior_var = inrange * 1/sum(inrange);

% endpoint distribution and expected gain function
% calculation of these functions are very costly
% these are invariant across subjects so taken out from Ln 79-82
eg1 = nan(length(Y),length(v_hut)); eg2 = nan(length(Y),length(v_hut));

tic
for iy = 1:length(Y)
    endpoint_dist_hut = exp(-(y - Y(iy) + 0).^2./(2*v_hut')) ./ (2*pi*v_hut').^(1/2);
    endpoint_dist_hut = endpoint_dist_hut ./ sum(endpoint_dist_hut,2);
    
    eg1(iy,:) = sum(endpoint_dist_hut .* gainFunc,2)'; % for all possible estimates of variance
    eg2(iy,:) = sum(endpoint_dist_hut .* gainFunc2,2)';
end
toc

save egfunc.mat  eg1 eg2 y v_hut Y prior_var -v7.3 -nocompression
