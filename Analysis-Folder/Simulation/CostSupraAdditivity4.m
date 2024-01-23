clear all; close all;

x = -20:0.05:25; %cm in steps of 0.5 mm
y = -15:0.05:15; 

[xx yy] = meshgrid(x,y);

s = mean([2.3 3.8]); % horizontal error s0 = 0.23-0.38 cm across subjects
Sigma = [s 0; 0 s]; %1sd~10mm=1cm
nrep = 100;
reward = 1;
superadd = @(p) 1*p + 0.063;
simresult = [];
xc = (max(x) + min(x))/2; yc = (max(y) + min(y))/2;
height = 50; % 5cm

% left target
width = Sigma(1)*1.5; % width of rectangles = s0, 1.5*s0, 2*s0
d = 0.6*width; % gap widths were 0.4, 0.6, 0.9 or 1.35 * the width of the corresponding rectangles.

T1(1) = xc-d/2-width ; T1(2) = xc-d/2 ; % left end right end
T1(3) = yc+height/2; T1(4) = yc-height/2; % top bottom

T2(1) = xc+d/2 ; T2(2) = xc+d/2+width ;
T2(3) = yc+height/2; T2(4) = yc-height/2;

Target1 = []; Target2 = []; 
Target1 = xx >= T1(1) & xx < T1(2) & yy >= T1(4) & yy <= T1(3);
Target2 = xx >= T2(1) & xx < T2(2) & yy >= T2(4) & yy <= T2(3);

% right target
width2 = Sigma(1)*0.84;  % adjusted by a stiar-case
T3(1) = xc-width2; T3(2) = xc+width2; 
T3(3) = yc+height/2; T3(4) = yc-height/2; % top bottom
Target3 = []; 
Target3 = xx >= T3(1) & xx < T3(2) & yy >= T3(4) & yy <= T3(3);

% endpoint
mux = xc; muy = yc;
samples = ndRandn([mux muy]', Sigma.^2, nrep)'; 
load('s3');
[s1, tmp]= sort(samples(:,1));
samples = [s1, samples(tmp,2)];
samples(50,1) = samples(50,1)+3;
samples(49,1) = samples(49,1)-1.5;


hitT1 = samples(:,1) >= T1(1) & samples(:,1) < T1(2) & samples(:,2) >= T1(4) & samples(:,2) <= T1(3);
hitT2 = samples(:,1) >= T2(1) & samples(:,1) < T2(2) & samples(:,2) >= T2(4) & samples(:,2) <= T2(3);
hitT3 = samples(:,1) >= T3(1) & samples(:,1) < T3(2) & samples(:,2) >= T3(4) & samples(:,2) <= T3(3);

% probability
probT1 = sum(hitT1) / nrep ;
probT2 = sum(hitT2) / nrep ;
prob = sum(hitT1 + hitT2) / nrep ;

probsub = superadd(prob);
result = [prob, probT1, probT2, probsub]  % P(L)=P(R)=0.3


probT3 = sum(hitT3) / nrep ;
logodds = @(p) log(p./(1-p));
pm = [0.88 0.76];
probT3_llo = pm(1) * logodds(probT3) + (1-pm(1)) * logodds(pm(2));
probT3_llo = exp(probT3_llo) ./ (1 + exp(probT3_llo)) ;
result2 = [probT3 probT3_llo]  % P(S)=0.6


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
% save s3 samples

figure;
rectangle('Position',[T3(1) T3(4) width2*2 height]); hold on 
plot(samples(:,1),samples(:,2),'k.', 'MarkerSize',11);
plot(mux, muy, 'gs', 'MarkerSize',6); 
axis('equal');
xlim([min(x) max(x)]); ylim([min(y) max(y)]);
box off;

pos(3) = 600; pos(4) = 600;
set(gcf, 'Position', pos);





