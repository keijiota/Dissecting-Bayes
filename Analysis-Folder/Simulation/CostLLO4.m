clear all; close all;

x = -7.5:0.01:7.5; %cm in steps of 0.5 mm
y = -5:0.01:5;
[xx yy] = meshgrid(x,y); 

Sigma = [1 0; 0 1]; %1sd~10mm=1cm
nrep = 1000;
reward = 1;
reward2 = 2;

logodds = @(p) log(p./(1-p));
pm = [0.88 0.72];

xc = (max(x) + min(x))/2; yc = (max(y) + min(y))/2;
height1 = 8;  height2 = 8;

% right target
Width3 = [0.5, 1, 1.5]*Sigma(1);

for ir = 1:3 
    width3 = Width3(ir); 
    simresult = []; simresult2 = [];

for is1 = 0:0.01:1
    % left target
    width1 = width3*is1;
    for is2 = 0:0.01:1
        width2 = width3*is2;
        
        if width1 + width2 <= width3
            
            T1(1) = xc-width1/2 ; T1(2) = xc+width1/2 ; % left end, right end
            T1(3) = yc+height1/2; T1(4) = yc-height1/2; % top, bottom
            
            T2L(1) = T1(1)-width2/2 ; T2L(2) = T1(1);
            T2L(3) = yc+height1/2;  T2L(4) = yc-height1/2;
            
            T2R(1) = T1(2) ; T2R(2) = T1(2)+width2/2;
            T2R(3) = yc+height1/2; T2R(4) = yc-height1/2;
            
            T3(1) = xc-width3; T3(2) = xc+width3 ;
            T3(3) = yc+height2/2; T3(4) = yc-height2/2;
            
            Target1 = []; Target2L = []; Target2R = []; Target3 = [];
%             for ix = 1:length(x)
%                 for iy = 1:length(y)
%                     if x(ix) >= T1(1) && x(ix) < T1(2) && y(iy) >= T1(4) && y(iy) <= T1(3)
%                         Target1(iy,ix) = 1;
%                     else
%                         Target1(iy,ix) = 0;
%                     end
%                     
%                     if x(ix) >= T2L(1) && x(ix) < T2L(2) && y(iy) >= T2L(4) && y(iy) <= T2L(3)
%                         Target2L(iy,ix) = 1;
%                     else
%                         Target2L(iy,ix) = 0;
%                     end
%                     
%                     if x(ix) > T2R(1) && x(ix) <= T2R(2) && y(iy) >= T2R(4) && y(iy) <= T2R(3)
%                         Target2R(iy,ix) = 1;
%                     else
%                         Target2R(iy,ix) = 0;
%                     end
%                     
%                     if x(ix) >= T3(1) && x(ix) <= T3(2) && y(iy) >= T3(4) && y(iy) <= T3(3)
%                         Target3(iy,ix) = 1;
%                     else
%                         Target3(iy,ix) = 0;
%                     end
%                     
%                 end
%             end

            Target1 = xx >= T1(1) & xx < T1(2) & yy >= T1(4) & yy <= T1(3); 
            Target2L = xx >= T2L(1) & xx < T2L(2) & yy >= T2L(4) & yy <= T2L(3);
            Target2R = xx > T2R(1) & xx <= T2R(2) & yy >= T2R(4) & yy <= T2R(3);
            Target3 = xx >= T3(1) & xx <= T3(2) & yy >= T3(4) & yy <= T3(3);
                       
            mux = xc; muy = yc;
            s = ndRandn([mux muy]', Sigma.^2, nrep)';
            
            hitT1 = s(:,1) >= T1(1) & s(:,1) < T1(2) & s(:,2) >= T1(4) & s(:,2) <= T1(3);
            hitT2L = s(:,1) >= T2L(1) & s(:,1) < T2L(2) & s(:,2) >= T2L(4) & s(:,2) <= T2L(3);
            hitT2R = s(:,1) > T2R(1) & s(:,1) <= T2R(2) & s(:,2) >= T2R(4) & s(:,2) <= T2R(3);
            hitT3 = s(:,1) >= T3(1) & s(:,1) <= T3(2) & s(:,2) >= T3(4) & s(:,2) <= T3(3);
            
            probT1 = sum(hitT1) / nrep ;
            probT2R = sum(hitT2L) / nrep ;
            probT2L = sum(hitT2R) / nrep ;
            probT2 = sum(hitT2R + hitT2L) / nrep ;
            probT3 = sum(hitT3) / nrep ;
            
            
            probT1_llo = pm(1) * logodds(probT1) + (1-pm(1)) * logodds(pm(2)); probT1_llo = exp(probT1_llo) ./ (1 + exp(probT1_llo));
            probT2_llo = pm(1) * logodds(probT2) + (1-pm(1)) * logodds(pm(2)); probT2_llo = exp(probT2_llo) ./ (1 + exp(probT2_llo));
            probT3_llo = pm(1) * logodds(probT3) + (1-pm(1)) * logodds(pm(2)); probT3_llo = exp(probT3_llo) ./ (1 + exp(probT3_llo));
            
            egLeft = reward2*probT1 + reward*probT2;
            egRight = reward*probT3;
            egLeft_llo = reward2*probT1_llo + reward*probT2_llo;
            egRight_llo = reward*probT3_llo;
            
            [eg chopt] = max([egLeft, egRight]);
            [eg_llo chact] = max([egLeft_llo, egRight_llo]);
            
            if chact == 1
                eg_act =  egLeft;
            else
                eg_act =  egRight;
            end
            if chopt == 1
                eg_opt =  egLeft;
            else
                eg_opt =  egRight;
            end
            
            [st tmp] = max(s(:,1)) ;                  
            chsalient = max(s(tmp,1)) <= T2R(2) & s(tmp,2) >= T2R(4) & s(tmp,2) <= T2R(3) ; 
            
            if chsalient == 1
                eg_salient =  egLeft;
            else
                eg_salient =  egRight;
            end            

            simresult = [simresult; eg_act / eg_opt, egLeft, egRight, egLeft_llo, egRight_llo, chact, chopt];
            simresult2 = [simresult2; eg_salient / eg_opt, egLeft, egRight, chsalient, chopt];
            
            
%             subplot(1,2,1)
%             TargetLeft = 2*Target1 + Target2L + Target2R;
%             mesh(x,y,TargetLeft); view([0 90]);
%             axis('equal');
%             xlim([min(x) max(x)]); ylim([min(y) max(y)]);
%             
%             subplot(1,2,2)            
%             TargetRight = Target3;
%             mesh(x,y,TargetRight); view([0 90]);
%             axis('equal');
%             xlim([min(x) max(x)]); ylim([min(y) max(y)]);
            
        end
    end
    
end

% pos(3) = 650; pos(4) = 600;
% set(gcf, 'Position', pos);

Simresult{1,ir} = simresult;
Simresult2{1,ir} = simresult2;

end

mean(Simresult{1,1})
mean(Simresult{1,2})
mean(Simresult{1,3})






