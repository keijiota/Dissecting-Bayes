clear ; close all ;

y = -300:1:500; v_hut = [0:0.001:0.1 0.2:0.1:1 2:1:5000]; m_hut = -10:0.1:10; Y = 1:10:180;
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
gainfunc = gainFunc; N1 = 30; N2 = 5;

%---- generative proccess -----
population_mean = [0 0];
population_var  = [10 20].^2;

% population pdf at the i-th trial, stimulus distribution
pop_mean = rand(1,1) * diff(population_mean) + population_mean(1);
pop_var = rand(1,1) * diff(population_var) + population_var(1);

% measurment distribution of population mean/sigma P(sample sigma | population sigma)
% mean_distribution = exp(-(m_hut-pop_mean).^2/(2*1^2)) / sqrt(2*pi*1^2); mean_distribution = mean_distribution/sum(mean_distribution);
var_distribution = (N1-1)/(2^((N1-1)/2)*gamma((N1-1)/2)*pop_var) .* ((v_hut*(N1-1))/pop_var).^((N1-1)/2-1) .* exp(-v_hut*(N1-1)/(2*pop_var));  var_distribution = var_distribution/sum(var_distribution);

% generate a sample
sample = randn(N1,1)* sqrt(pop_var) + pop_mean;

%---- inference proccess -----
% mearuement
sample_mean = mean(sample);
sample_var1 = var(sample,0); % sum((sample-mean(sample)).^2)/(N-1)
sample_var1 = 200;
sample_var2 = 200;

% prior probability
inrange =  m_hut >= population_mean(1) & m_hut <= population_mean(end) ;
prior_mean = inrange * 1/sum(inrange);

inrange =  v_hut >= population_var(1) & v_hut <= population_var(end) ;
prior_var = inrange * 1/sum(inrange);

% likelihood of population mean/sigma P(sample sigma | population sigma)
likelihood_var1  = (N1-1)/(2^((N1-1)/2)*gamma((N1-1)/2)*sample_var1) .* ((v_hut*(N1-1))/sample_var1).^((N1-1)/2-1) .* exp(-v_hut*(N1-1)/(2*sample_var1));  likelihood_var1 = likelihood_var1/sum(likelihood_var1);
likelihood_var2  = (N2-1)/(2^((N2-1)/2)*gamma((N2-1)/2)*sample_var2) .* ((v_hut*(N2-1))/sample_var2).^((N2-1)/2-1) .* exp(-v_hut*(N2-1)/(2*sample_var2));  likelihood_var2 = likelihood_var2/sum(likelihood_var2);

% posterior of population mean/sigma P(population sigma | sample sigma)
posterior_var1 = prior_var .* likelihood_var1; posterior_var1 = posterior_var1/sum(posterior_var1);
posterior_var2 = prior_var .* likelihood_var2; posterior_var2 = posterior_var2/sum(posterior_var2);
sum(posterior_var1.*v_hut)
sum(posterior_var2.*v_hut)


vmax = 600;
subplot(3,4,1); plot(v_hut,prior_var,'r');
myfig; xticks([100 400 vmax]); xlim([0 vmax]); xlabel('Estimates of sample variance');
subplot(3,4,5); plot(v_hut,likelihood_var1,'r'); hold on
myfig; xticks([100 400 vmax]); xlim([0 vmax]); xlabel('Estimates of sample variance');
subplot(3,4,5); plot(v_hut,likelihood_var2,'b'); hold on
lineplot(sample_var1, 'v','r-'); lineplot(sample_var2, 'v','b-');
subplot(3,4,9); plot(v_hut,posterior_var1,'r'); hold on
subplot(3,4,9); plot(v_hut,posterior_var2,'b');
myfig; xticks([100 400 vmax]); xlim([0 vmax]); xlabel('Estimates of sample variance');

%---- calculation of expected gain -----
Y = 150;
endpoint_dist_hut = exp(-(y - Y + 0).^2./(2*v_hut')) ./ sqrt(2*pi*v_hut');
endpoint_dist_hut = endpoint_dist_hut ./ sum(endpoint_dist_hut,2);

tmp1 = 50;  tmp1 = find(abs(v_hut - tmp1) < 0.01);
tmp2 = 250; tmp2 = find(abs(v_hut - tmp2) < 0.01);
tmp3 = 550; tmp3 = find(abs(v_hut - tmp3) < 0.01);
ymax = 230; a = max(endpoint_dist_hut(tmp1,:)); b = 0.001;
subplot(3,4,2); xticks(0:100:200); plot(y, endpoint_dist_hut(tmp1,:),'r'); hold on
myfig; xticks([90 180 ymax]); xlim([0 ymax]); ylim([0 a+b]);  xlabel('Endpoint');
lineplot(Y, 'v','r-');
subplot(3,4,6); plot(y, endpoint_dist_hut(tmp2,:),'r'); hold on
myfig; xticks([90 180 ymax]); xlim([0 ymax]); ylim([0 a+b]); xlabel('Endpoint');
lineplot(Y, 'v','r-');
subplot(3,4,10); plot(y, endpoint_dist_hut(tmp3,:),'r'); hold on
myfig; xticks([90 180 ymax]); xlim([0 ymax]); ylim([0 a+b]); xlabel('Endpoint');
lineplot(Y, 'v','r-');

egfunction = endpoint_dist_hut .* gainfunc;

a = max(egfunction(tmp1,:)); b = 0.01;
subplot(3,4,4); plot(y,gainfunc,'r');
myfig; xticks([90 180 ymax]); xlim([0 ymax]); xlabel('Endpoint');
subplot(3,4,3); plot(y, egfunction(tmp1,:),'r');
myfig; xticks([90 180 ymax]); xlim([0 ymax]); ylim([0 a+b]); xlabel('Endpoint');
subplot(3,4,7); plot(y, egfunction(tmp2,:),'r');
myfig; xticks([90 180 ymax]); xlim([0 ymax]); ylim([0 a+b]); xlabel('Endpoint');
subplot(3,4,11); plot(y, egfunction(tmp3,:),'r');
myfig; xticks([90 180 ymax]); xlim([0 ymax]); ylim([0 a+b]); xlabel('Endpoint');

eg = sum(endpoint_dist_hut .* gainfunc,2)'; % for all possible estimates of variance, integral for endpoint
subplot(3,4,8); plot(v_hut, eg,'r');
myfig; xticks([100 400 vmax]); xlim([0 vmax]); ylim([60 100]); xlabel('Estimates of sample variance');

expectedeg = eg.*posterior_var1; % scaled by a posterior of population variance
expectedeg2 = eg.*posterior_var2;
subplot(3,4,12); plot(v_hut, expectedeg,'r'); hold on
subplot(3,4,12); plot(v_hut, expectedeg2,'b');
myfig; xticks([100 400 vmax]); xlim([0 vmax]); xlabel('Estimates of sample variance');

ExpectedExpectedGain = nansum(eg.*posterior_var1) % expected expected gain taking a posterior into acount, integral for estimates of variance
ExpectedExpectedGain = nansum(eg.*posterior_var2)

%
% figure
% plot(y, endpoint_dist_hut(2,:))
% plot(y, endpoint_dist_hut(2,:) .* gainfunc)
sum(endpoint_dist_hut(2,:) .* gainfunc)

set(gcf, 'Position', [0 0 1300 470]);


% likelihood of population mean/sigma P(sample sigma | population sigma)
sample_var1 = 100; sample_var2 = 100;
likelihood_var1  = (N1-1)/(2^((N1-1)/2)*gamma((N1-1)/2)*sample_var1) .* ((v_hut*(N1-1))/sample_var1).^((N1-1)/2-1) .* exp(-v_hut*(N1-1)/(2*sample_var1));  likelihood_var1 = likelihood_var1/sum(likelihood_var1);
likelihood_var2  = (N2-1)/(2^((N2-1)/2)*gamma((N2-1)/2)*sample_var2) .* ((v_hut*(N2-1))/sample_var2).^((N2-1)/2-1) .* exp(-v_hut*(N2-1)/(2*sample_var2));  likelihood_var2 = likelihood_var2/sum(likelihood_var2);

sample_var3 = 300; sample_var4 = 300;
likelihood_var3  = (N1-1)/(2^((N1-1)/2)*gamma((N1-1)/2)*sample_var3) .* ((v_hut*(N1-1))/sample_var3).^((N1-1)/2-1) .* exp(-v_hut*(N1-1)/(2*sample_var3));  likelihood_var3 = likelihood_var3/sum(likelihood_var3);
likelihood_var4  = (N2-1)/(2^((N2-1)/2)*gamma((N2-1)/2)*sample_var4) .* ((v_hut*(N2-1))/sample_var4).^((N2-1)/2-1) .* exp(-v_hut*(N2-1)/(2*sample_var4));  likelihood_var4 = likelihood_var4/sum(likelihood_var4);

figure
for sample_var1 = [50 150 250 350 450]
    sample_var2 = sample_var1;
    likelihood_var1  = (N1-1)/(2^((N1-1)/2)*gamma((N1-1)/2)*sample_var1) .* ((v_hut*(N1-1))/sample_var1).^((N1-1)/2-1) .* exp(-v_hut*(N1-1)/(2*sample_var1));  likelihood_var1 = likelihood_var1/sum(likelihood_var1);
    likelihood_var2  = (N2-1)/(2^((N2-1)/2)*gamma((N2-1)/2)*sample_var2) .* ((v_hut*(N2-1))/sample_var2).^((N2-1)/2-1) .* exp(-v_hut*(N2-1)/(2*sample_var2));  likelihood_var2 = likelihood_var2/sum(likelihood_var2);
    subplot(1,2,1)
    plot(v_hut,likelihood_var1); hold on
    subplot(1,2,2)
    plot(v_hut,likelihood_var2); hold on
end
vmax = 1000;
subplot(1,2,1)
xticks([100 400 vmax]); xlim([0 vmax]); xlabel('Estimates of sample variance');
ylim([0 0.035]); yticks(0:0.01:0.03); myfig;

subplot(1,2,2)
xticks([100 400 vmax]); xlim([0 vmax]); xlabel('Estimates of sample variance');
ylim([0 0.035]); yticks(0:0.01:0.03);  myfig;

