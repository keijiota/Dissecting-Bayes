clear ; close all;

rng('shuffle');

load sim_weight;
load B ;
n = 100;
for con = 1:4
    
    b = meanW{1,con}' ;
    y = s(1:n,con);
    X = xx{1,con};
    X = X(1:n,:);
    Xc = [X, ones(length(y),1)];
    N = size(X,2);
    % [~, vif] = multicollinearity(X)
    
    b1 = Xc \ y; b1 = b1(1:end-1);
    b2 = inv(Xc'*Xc)*Xc'*y; b2 = b2(1:end-1);
    det(X'*X)
    
    mdl = fitglm(X,y);
    b3 = mdl.Coefficients.Estimate;
    b3 = b3(2:end);
    
    [b4,fitinfo4]  = lasso(X,y,'CV',10,'Standardize',true);
%     figure; lassoPlot(b4,fitinfo4,'PlotType','Lambda','XScale','log') ;
    figure; lassoPlot(b4,fitinfo4,'PlotType','CV') ;
    tmp = fitinfo4.IndexMinMSE ;  b4a = b4(:,tmp) ;
    tmp = fitinfo4.Index1SE ;     b4b = b4(:,tmp) ;
    
    [b5,fitinfo]  = lasso(X,y,'CV',10,'Standardize',false, 'Alpha',0.01);
    tmp = fitinfo.IndexMinMSE ;  b5 = b5(:,tmp) ;
    
    k = 10; standardize = 0; figureon = 1;
    [Bridge, fitinfo] = myridge(y,X,k,standardize,figureon);
    b6a = Bridge(:,fitinfo.indexminmse);
    b6b = Bridge(:,fitinfo.index1se);
    
    I = eye(size(Xc,2));
    b7 = inv(Xc'*Xc + fitinfo.lambdaMin*I)*Xc'*y;  b7 = b7(1:end-1) ;
    
    B{1,con} = [b,b1,b2,b3,b4a,b4b,b5,b6a,b6b,b7];
    
    figure(30);
    subplot(2,2,con); hold on
    plot(1:N, b1,'r-'); 
    plot(1:N, b4a,'g-'); plot(1:N, b6a,'b-');
    plot(1:N, b4b,'g--'); plot(1:N, b6b,'b--');
    subplot(2,2,con); plot(1:N, b,'k-'); ylim([min(b)-0.05, max(b)+0.05]);
    
end



