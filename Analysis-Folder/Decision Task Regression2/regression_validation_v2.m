clear ; close all;

load sim_weight;
load B ;
n = 100;
for con = 1
    
    b = meanW{1,con}' ;
    y = s(1:n,con);
    X = xx{1,con};
    X = X(1:n,:);
    Xc = [X, ones(length(y),1)];
    N = size(X,2);
    
    for irep = 1:5
        [bl,fitinfol]  = lasso(X,y,'CV',10,'Standardize',false,'Alpha',1);
%             figure; lassoPlot(bl,fitinfol,'PlotType','Lambda','XScale','log') ;
%              figure; lassoPlot(bl,fitinfol,'PlotType','CV') ;
             
        minmse(irep,1) = min(fitinfol.MSE);
        lambdaminmse(irep,1) = fitinfol.LambdaMinMSE;
        lambdaonese(irep,1) = fitinfol.Lambda1SE;
        blasso_minmse(irep,:) = bl(:,fitinfol.IndexMinMSE);
        blasso_onese(irep,:) = bl(:,fitinfol.Index1SE);
        
        k = 10; standardize = 0; figureon = 0;
        [Br, fitinfor] = myridge(y,X,k,standardize,figureon);
        minmse(irep,2) = min(fitinfor.MSE);
        lambdaminmse(irep,2) = fitinfor.LambdaMinMSE;
        lambdaonese(irep,2) = fitinfor.Lambda1SE;
        bridge_minmse(irep,:) = Br(:,fitinfor.IndexMinMSE);
        bridge_onese(irep,:) = Br(:,fitinfor.Index1SE);
    end
    
    figure;
    plot(1:N,b,'k-','linewidth',1.5); hold on
    plot(1:N,bridge_minmse,'b:'); plot(1:N,mean(bridge_minmse),'b-','linewidth',2);
    plot(1:N,blasso_minmse,'g:'); plot(1:N,mean(blasso_minmse),'g-','linewidth',2);
    ylim([min(b)-0.05 max(b)+0.05]);
    
    figure;
    plot(1:N,b,'k-','linewidth',1.5); hold on
    plot(1:N,bridge_onese,'b:'); plot(1:N,mean(bridge_onese),'b-','linewidth',2);
    plot(1:N,blasso_onese,'g:'); plot(1:N,mean(blasso_onese),'g-','linewidth',2);
    ylim([min(b)-0.05 max(b)+0.05]);
        
end






