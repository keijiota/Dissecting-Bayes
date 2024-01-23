clear all; close all;

x = -20:0.05:25; %mm in steps of 0.05 mm
y = -15:0.05:15; 
[xx yy] = meshgrid(x,y);

s  = mean([2.3 3.8]); % horizontal error s0 = 0.23-0.38 cm across subjects
Sigma = [s 0; 0 s]; %1sd~10mm=1cm
xc = (max(x) + min(x))/2; yc = (max(y) + min(y))/2;
height = 50; % 5cm

% left target
width = Sigma(1)*1.5; % width of rectangles = s0, 1.5*s0, 2*s0
d = 0.75*width; % gap widths were 0.4, 0.6, 0.9 or 1.35 * the width of the corresponding rectangles.

T1(1) = xc-d/2-width ; T1(2) = xc-d/2 ; % left end right end
T1(3) = yc+height/2; T1(4) = yc-height/2; % top bottom

T2(1) = xc+d/2 ; T2(2) = xc+d/2+width ;
T2(3) = yc+height/2; T2(4) = yc-height/2;

Target1 = []; Target2 = []; 
Target1 = xx >= T1(1) & xx < T1(2) & yy >= T1(4) & yy <= T1(3);
Target2 = xx >= T2(1) & xx < T2(2) & yy >= T2(4) & yy <= T2(3);

% right target
width2 = Sigma(1)*0.73;  % adjusted to double
T3(1) = xc-width2; T3(2) = xc+width2; 
T3(3) = yc+height/2; T3(4) = yc-height/2; % top bottom
Target3 = []; Target3 = xx >= T3(1) & xx < T3(2) & yy >= T3(4) & yy <= T3(3);

% probability from pdf (~= many emprical samples)
mux = xc; muy = yc; nrep = 100000000;
samples = ndRandn([mux muy]', Sigma.^2, nrep)'; 

hitT1 = samples(:,1) >= T1(1) & samples(:,1) < T1(2) & samples(:,2) >= T1(4) & samples(:,2) <= T1(3);
hitT2 = samples(:,1) >= T2(1) & samples(:,1) < T2(2) & samples(:,2) >= T2(4) & samples(:,2) <= T2(3);
probT1 = sum(hitT1) / nrep ;
probT2 = sum(hitT2) / nrep ;
probDouble   = sum(hitT1 + hitT2) / nrep ;
superadd = @(p) 1*p + 0.063;
probDoublesub = superadd(probDouble);
[probDouble, probT1, probT2, probDoublesub] % P(L)=P(R)=0.2673

hitT3 = samples(:,1) >= T3(1) & samples(:,1) < T3(2) & samples(:,2) >= T3(4) & samples(:,2) <= T3(3);
probSingle = sum(hitT3) / nrep ;
logodds = @(p) log(p./(1-p));
pm = [0.88 0.76];
probSingle_sub = pm(1) * logodds(probSingle) + (1-pm(1)) * logodds(pm(2));
probSingle_sub = exp(probSingle_sub) ./ (1 + exp(probSingle_sub)) ;
[probSingle probSingle_sub]  % P(S)=0.5346


% right target adjusted to double
width2adj = Sigma(1)*0.8373;  % adjusted by a stiar-case
T3adj(1) = xc-width2adj; T3adj(2) = xc+width2adj; 
T3adj(3) = yc+height/2; T3adj(4) = yc-height/2; % top bottom
Target3adj = []; Target3adj = xx >= T3adj(1) & xx < T3adj(2) & yy >= T3adj(4) & yy <= T3adj(3);
hitT3adj = samples(:,1) >= T3adj(1) & samples(:,1) < T3adj(2) & samples(:,2) >= T3adj(4) & samples(:,2) <= T3adj(3);
probSingleAdj = sum(hitT3adj) / nrep;
[probSingleAdj width2*2 width2adj*2 width2adj*2/(width2*2)] %P(S_adj)=0.5976


% endpoint
load('s3');

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

figure;
rectangle('Position',[T3(1) T3(4) width2*2 height]); hold on 
rectangle('Position',[T3adj(1) T3adj(4) width2adj*2 height]);
plot(samples(:,1),samples(:,2),'k.', 'MarkerSize',11);
plot(mux, muy, 'gs', 'MarkerSize',6); 
axis('equal');
xlim([min(x) max(x)]); ylim([min(y) max(y)]);
box off;
pos(3) = 600; pos(4) = 600;
set(gcf, 'Position', pos);

% transformation
figure;  hold on
x = 0:0.01:1;
y = superadd(x);
plot(probDouble,probDoublesub, 'ko', 'MarkerFaceColor','w', 'MarkerSize',7) ; hold on
plot(0:0.1:1.2,0:0.1:1.2,'k--');
plot(x, y, 'k-');
axis('square');
xlim([0 1.05]); ylim([0 1.05]);
xticks(0:0.2:1); yticks(0:0.2:1);
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);
box off; 
pos(3) = 200; pos(4) = 200;
set(gcf, 'Position', pos);

figure; hold on
y = pm(1) * logodds(x) + (1-pm(1)) * logodds(pm(2));
y = exp(y) ./ (1 + exp(y)) ;
plot(probSingle,probSingle_sub, 'ko', 'MarkerFaceColor','w', 'MarkerSize',7) ; hold on
plot(0:0.1:1.2,0:0.1:1.2,'k--');
plot(x, y, 'k-');
axis('square');
xlim([0 1.05]); ylim([0 1.05]);
xticks(0:0.2:1); yticks(0:0.2:1);
set(gca, 'Fontname', 'Arial Regular', 'Fontsize', 12, 'linewidth', 1,'TickLength',[0.025 0]);
box off; 
pos(3) = 200; pos(4) = 200;
set(gcf, 'Position', pos);










