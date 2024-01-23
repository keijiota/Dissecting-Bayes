clear ; close all ;

load('dat_decision');
load suffdata
% figure; hold on 
tmp = [1 3 2 4]; % [30,0], [5,0], [30,-500], [5,-500]
for con = 1:4
    % 17 subs * 50 trials * 4 conditions
    optaim = OptAim_Sufficient_G(:,:,con);
    samplevar  = SampleNoiseY(:,:,con).^2;
    optaim = reshape(optaim, 17*50,1);
    samplevar = reshape(samplevar, 17*50,1);
    [samplevar idx]= sort(samplevar);
    optaim = optaim(idx);  
    OPTAim(:,con) = optaim;
    SampleVar(:,con) = samplevar;
end

tmp = [1 3 2 4]; % [30,0], [5,0], [30,-500], [5,-500]
OPTAim = OPTAim(:,tmp);
SampleVar = SampleVar(:,tmp);

ms = 5.5; cs = 10; ms2 = 10;
figure; hold on 
C1 = [204 247 255]./255 ; C11 = [80 180 255]./255 ; 
C3 = [179 179 255]./255; C33 = [25 25 255]./255;
C2 = [255 221 200]./255; C22 = [255 140 0]./255;
C4 = [255 200 200]./255; C44 = [200 0 0]./255; 

plot(SampleVar(:,1), OPTAim(:,1), '-', 'Color', C11,'LineWidth',2); hold on
plot(SampleVar(:,2), OPTAim(:,2), '-', 'Color', C33,'LineWidth',2); hold on
plot(SampleVar(:,3), OPTAim(:,3), '-', 'Color', C22,'LineWidth',2); hold on
plot(SampleVar(:,4), OPTAim(:,4), '-', 'Color', C44,'LineWidth',2); hold on
yline(180,'g-','LineWidth',1.5)
yticks(120:10:180); ylim([120 185]); 
xlim([0 1800]); xticks(0:250:1800);
myfigAI2(12,0.025, 500,400)

