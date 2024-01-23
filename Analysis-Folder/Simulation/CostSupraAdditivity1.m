clear all; close all;

x = -7.5:0.05:7.5; %cm in steps of 0.5 mm
y = -5:0.05:5;

Sigma = [1 0; 0 1]; %1sd~10mm=1cm
nrep = 1000000; 
reward = 1;
superadd = @(p) 1*p + 0.063; 
simresult = [];
xc = (max(x) + min(x))/2; yc = (max(y) + min(y))/2;
height1 = 8;  height2 = 8;

for is = 1:3
    if is == 1
        % right target
        width2 = Sigma(1)*1; d2 = 0; 
        % left target
        width1 = width2*0.47; d1 = 0.25;        
    elseif is == 2
        width2 = Sigma(1)*1.5; d2 = 0;        
        width1 = width2*0.51; d1 = 0.25;        
    elseif is == 3
        width2 = Sigma(1)*2.5; d2 = 0;                
        width1 = width2*0.65; d1 = 0.25;        
    end

T1(1) = xc-d1-width1 ; % left end
T1(2) = xc-d1 ; % right end
T1(3) = yc+height1/2; % top
T1(4) = yc-height1/2; % bottom

T2(1) = xc+d1 ; T2(2) = xc+d1+width1 ;
T2(3) = yc+height1/2; T2(4) = yc-height1/2;

T3(1) = xc+d2 - width2/2; T3(2) = xc+d2+width2/2 ;
T3(3) = yc+height2/2; T3(4) = yc-height2/2;

Target1 = []; Target2 = []; Target3 = [];
for ix = 1:length(x)
    for iy = 1:length(y)
        if x(ix) >= T1(1) && x(ix) <= T1(2) && y(iy) >= T1(4) && y(iy) <= T1(3)
            Target1(iy,ix) = 1;
        else
            Target1(iy,ix) = 0;
        end
        
        if x(ix) >= T2(1) && x(ix) <= T2(2) && y(iy) >= T2(4) && y(iy) <= T2(3)
            Target2(iy,ix) = 1;
        else
            Target2(iy,ix) = 0;
        end
        
        if x(ix) >= T3(1) && x(ix) <= T3(2) && y(iy) >= T3(4) && y(iy) <= T3(3)
            Target3(iy,ix) = 1;
        else
            Target3(iy,ix) = 0;
        end
        
    end
end


mux = xc; muy = yc;
s = ndRandn([mux muy]', Sigma.^2, nrep)';

hitT1 = s(:,1) >= T1(1) & s(:,1) <= T1(2) & s(:,2) >= T1(4) & s(:,2) <= T1(3);
hitT2 = s(:,1) >= T2(1) & s(:,1) <= T2(2) & s(:,2) >= T2(4) & s(:,2) <= T2(3);
hitT3 = s(:,1) >= T3(1) & s(:,1) <= T3(2) & s(:,2) >= T3(4) & s(:,2) <= T3(3);

probT1 = sum(hitT1) / nrep ;
probT2 = sum(hitT2) / nrep ;
probT3 = sum(hitT3) / nrep ;

probLeft = sum(hitT1 + hitT2) / nrep ;
probRight = sum(hitT3) / nrep ;
probLeftsub = superadd(probLeft);

[eg chopt] = max([reward*probLeft, reward*probRight]); 
[eg_superadd chact] = max([reward*superadd(probLeft), reward*probRight]);

if chact == 1
   eg_act =  reward*probLeft; 
else 
   eg_act =  reward*probRight; 
end
if chopt == 1
   eg_opt =  reward*probLeft; 
else 
   eg_opt =  reward*probRight; 
end

simresult = [simresult; eg_act / eg_opt, probLeft, probRight, probLeftsub, chact,chopt];


subplot(3,2,is*2-1)
TargetLeft = Target1 + Target2;
mesh(x,y,TargetLeft); view([0 90]);
axis('equal');
xlim([min(x) max(x)]); ylim([min(y) max(y)]);

subplot(3,2,is*2)
TargetRight = Target3;
mesh(x,y,TargetRight); view([0 90]);
axis('equal');
xlim([min(x) max(x)]); ylim([min(y) max(y)]);

end
pos(3) = 650; pos(4) = 600;
set(gcf, 'Position', pos);

simresult










