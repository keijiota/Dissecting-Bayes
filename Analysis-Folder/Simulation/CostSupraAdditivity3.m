clear all; close all;

x = -22.5:0.05:22.5; %cm in steps of 0.5 mm
y = -15:0.05:15; 

[xx yy] = meshgrid(x,y);

s = mean([12, 23, 11, 20, 11].^(1/2));
Sigma = [s 0; 0 s]; %1sd~10mm=1cm
nrep = 100;
reward = 1;
superadd = @(p) 1*p + 0.063;
simresult = [];
xc = (max(x) + min(x))/2; yc = (max(y) + min(y))/2;
height = 11*2;

% left target
width = Sigma(1)*1;
d = 0.2*width;

T1(1) = xc-d/2-width ; T1(2) = xc-d/2 ; % left end right end
T1(3) = yc+height/2; T1(4) = yc-height/2; % top bottom

T2(1) = xc+d/2 ; T2(2) = xc+d/2+width ;
T2(3) = yc+height/2; T2(4) = yc-height/2;

Target1 = []; Target2 = []; 
Target1 = xx >= T1(1) & xx < T1(2) & yy >= T1(4) & yy <= T1(3);
Target2 = xx >= T2(1) & xx < T2(2) & yy >= T2(4) & yy <= T2(3);

% right target
width2 = Sigma(1)*1.01
T3(1) = xc-width2; T3(2) = xc+width2; 
T3(3) = yc+height/2; T3(4) = yc-height/2; % top bottom
Target3 = []; 
Target3 = xx >= T3(1) & xx < T3(2) & yy >= T3(4) & yy <= T3(3);

% endpoint
mux = xc; muy = yc;
samples = ndRandn([mux muy]', Sigma.^2, nrep)';
load('s2');

hitT1 = samples(:,1) >= T1(1) & samples(:,1) < T1(2) & samples(:,2) >= T1(4) & samples(:,2) <= T1(3);
hitT2 = samples(:,1) >= T2(1) & samples(:,1) < T2(2) & samples(:,2) >= T2(4) & samples(:,2) <= T2(3);
hitT3 = samples(:,1) >= T3(1) & samples(:,1) < T3(2) & samples(:,2) >= T3(4) & samples(:,2) <= T3(3);

% probability
probT1 = sum(hitT1) / nrep ;
probT2 = sum(hitT2) / nrep ;
prob = sum(hitT1 + hitT2) / nrep ;

probsub = superadd(prob);
result = [prob, probT1, probT2, probsub]


probT3 = sum(hitT3) / nrep ;
logodds = @(p) log(p./(1-p));
pm = [0.88 0.76];
probT3_llo = pm(1) * logodds(probT3) + (1-pm(1)) * logodds(pm(2));
probT3_llo = exp(probT3_llo) ./ (1 + exp(probT3_llo)) ;
result2 = [probT3 probT3_llo]

% TargetLeft = Target1 + Target2;
% mesh(x,y,TargetLeft); view([0 90]);

figure;
rectangle('Position',[T1(1) T1(4) width height]); hold on 
rectangle('Position',[T2(1) T2(4) width height]); 
plot(samples(:,1),samples(:,2),'k.', 'MarkerSize',11);
plot(mux, muy, 'gs', 'MarkerSize',6); 
axis('equal');
xlim([min(x) max(x)]); ylim([min(y) max(y)]);
box off;

pos(3) = 600; pos(4) = 600;
set(gcf, 'Position', pos);
% save s2 samples

figure;
rectangle('Position',[T3(1) T3(4) width2*2 height]); hold on 
plot(samples(:,1),samples(:,2),'k.', 'MarkerSize',11);
plot(mux, muy, 'gs', 'MarkerSize',6); 
axis('equal');
xlim([min(x) max(x)]); ylim([min(y) max(y)]);
box off;

pos(3) = 600; pos(4) = 600;
set(gcf, 'Position', pos);





