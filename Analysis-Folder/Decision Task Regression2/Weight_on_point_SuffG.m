clear ; close all;

load egfunc; 

rng('shuffle');
nrep = 1000;
s = nan(nrep,4); 

tic
for con = 1:4
    if con == 1
        N = 30; eg = eg1;
    elseif con == 2
        N = 30; eg = eg2;
    elseif con == 3
        N = 5; eg = eg1;
    else
        N = 5; eg = eg2;
    end
    
    W = nan(nrep,N); x = [];   
    parfor irep = 1:nrep
        sigma = rand*10+10;
        X = randn(1,N)*sigma;
        X = sort(X,2) ;
        x(irep,:) = X;
        
        samplesigma = std(X,[],2);
        [S] = cal_Sstar(samplesigma,N,eg,prior_var,v_hut,Y);
        s(irep,con) = S;
        
        dXi = 1; shiftS = nan(1,N);
        for ix = 1:N
            XX = X; XX(1,ix) = X(1,ix) + dXi ;
            shiftsamplesigma = std(XX,[],2) ;
            [S] = cal_Sstar(shiftsamplesigma,N,eg,prior_var,v_hut,Y);
            shiftS(1,ix) = S ;
        end
        
        dS_dXi = shiftS - s(irep,con) ; % dS/dXi = dXi*Wi => Wi = dS/dXi / dXi
        W(irep,:) = dS_dXi / dXi ;
    end
    xx{1,con} = x;
    % subplot(2,2,con); plot(1:N, mean(W),'k-');
    meanW{1,con} = mean(W);
end
toc

save sim_weight2 meanW xx s


% %%
% load sim_weight;
% % for con = 1:4
% % if con == 1
% %     N = 30; 
% % elseif con == 2
% %     N = 30; 
% % elseif con == 3
% %     N = 5;
% % else
% %     N = 5;
% % end    
% % b = meanW{1,con}' ;
% % figure(1); 
% % subplot(2,2,con); hold on
% % subplot(2,2,con); plot(1:N, b,'k-');
% % end
% 
% nrep = 10000; s = nan(nrep,4); xx = [];
% tic
% for con = 1:4
% if con == 1
%     N = 30; eg = eg1;
% elseif con == 2
%     N = 30; eg = eg2;
% elseif con == 3
%     N = 5; eg = eg1;
% else
%     N = 5; eg = eg2;
% end
% 
% x = [];  
% parfor irep = 1:nrep
%     sigma = rand*10+10;
%     X = randn(1,N)*sigma;
%     X = sort(X,2) ;
%     x(irep,:) = X;
%     
%     samplesigma = std(X,[],2);
%     [S] = cal_Sstar(samplesigma,N,eg,prior_var,v_hut,Y);
%     s(irep,con) = S;    
% end
% 
% xx{1,con} = x; 
% y = s(:,con); b = meanW{1,con}' ;
% X = xx{1,con}; Xc = [X, ones(nrep,1)];
% % [~, vif] = multicollinearity(X)
% 
% b1 = Xc \ y; b1 = b1(1:end-1);
% b2 = inv(Xc'*Xc)*Xc'*y; b2 = b2(1:end-1);
% det(X'*X)
% 
% mdl = fitglm(X,y);
% b3 = mdl.Coefficients.Estimate;
% b3 = b3(2:end);
% 
% [b4,fitinfo]  = lasso(X,y,'CV',10,'Standardize',true);
% tmp = fitinfo.IndexMinMSE ;  b4 = b4(:,tmp) ;
% 
% [b5,fitinfo]  = lasso(X,y,'CV',10,'Standardize',false, 'Alpha',0.01);
% tmp = fitinfo.IndexMinMSE ;  b5 = b5(:,tmp) ;
% 
% 
% k = 10; lamda = -5:1:8;
% mse_ik =[];
% indx = randi(nrep,nrep,k);
% for ik = 1:k
%     indx_ik = indx(:,ik);
%     X_train  = X(indx_ik(k+1:end),:);  X_validate  = X(indx_ik(1:k),:);
%     Xc_train = Xc(indx_ik(k+1:end),:); Xc_validate = Xc(indx_ik(1:k),:);
%     y_train  = y(indx_ik(k+1:end),1);  y_validate  = y(indx_ik(1:k),1);
%     
%     br = ridge(y_train,X_train,10.^lamda,0);
%     br = [br(2:end,:); br(1,:)];
%     y_hut = Xc_validate * br;
%     mse_ik(ik,:) = mean((y_validate - y_hut).^2);
% end
% mse = mean(mse_ik); se = std(mse_ik) ./ sqrt(k);
% 
% figure(2);
% subplot(1,2,1); loglog(10.^lamda, mse);
% subplot(1,2,2); semilogx(10.^lamda, br(1:end-1,:));
% 
% [~, tmp] = min(mse); minlamda = 10.^lamda(tmp);
% 
% b6 = ridge(y,X,minlamda,0); b6 = b6(2:end);
% 
% I = eye(size(Xc,2));
% b7 = inv(Xc'*Xc + minlamda*I)*Xc'*y;  b7 = b7(1:end-1);
% 
% B{1,con} = [b,b1,b2,b3,b4,b5,b6,b7];
% 
% figure(1); 
% subplot(2,2,con); hold on
% plot(1:N, b1,'r-');
% plot(1:N, b4,'g-'); plot(1:N, b6,'b-');
% subplot(2,2,con); plot(1:N, b,'k-');
% 
% end
% toc
% 
% save B B s xx
% 
% 
% % save B B
% 
% 
% 

