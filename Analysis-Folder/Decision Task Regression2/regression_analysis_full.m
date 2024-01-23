clear; close all ;

load('dat_decision');

for con = 1:4
    endpointorg = EndpointorgY(:,:,con);
    samplesY = SamplesY(:,:,con);
    orgnoise = OrgNoiseY(:,:,con);
    samplenoise = SampleNoiseY(:,:,con);
    optAim = OptAim_Sufficient_G(:,:,con);
    
    ssy = []; sy = []; ssamples = [];
    
    [N T] = size(endpointorg);
    for subi = 1:N
        for itrial = 1:T
            % each samples sorted in a order
            sy = samplesY(subi, itrial*31-30:itrial*31-1) ; %X1-X30
            if con == 3 || con == 4
                sy(6) = nan; % sixth point is feedback point
            end
            ssy = sort(sy);
            ssamples(itrial, :, subi) = ssy;
        end
    end
    
    for subi = 1:17
        endpoint = endpointorg(subi, :)';
        optaim = optAim(subi,:)';        
        sample = ssamples(:,:,subi);
        PopSD = orgnoise(subi,:);
        
        if con == 1 || con == 2
            orgsample = sample - endpoint ;
        else
            orgsample = sample(:,1:5) - endpoint ;
        end
        
        X = [orgsample, mean(orgsample')', std(orgsample')', PopSD'];
        % sorted raw data (not normalized), sample mean, sample SD, population SD
        % Variance is squared and not good for a linear model
        X = [orgsample];
        y = endpoint;
        y = optaim ; 
        nfeatures = size(X,2) ;        
        
        for ir = 1:4
            if ir == 1
                [b,fitinfo] = lasso(X,y,'CV',10,'Standardize',false,'Alpha',1);
%                 lassoPlot(b,fitinfo,'PlotType','Lambda','XScale','log') ;
%                 lassoPlot(b,fitinfo,'PlotType','CV') ;                              
            elseif ir == 2
                [b,fitinfo] = lasso(X,y,'CV',10,'Standardize',false,'Alpha',0.01);
%                 lassoPlot(b,fitinfo,'PlotType','Lambda','XScale','log') ;              
%                 lassoPlot(b,fitinfo,'PlotType','CV') ;                                              
            elseif ir == 3
                k = 10; standardize = 0; figureon = 0;
                [b, fitinfo] = myridge(y,X,k,standardize,figureon);
            else
                mdl = fitglm(X,y);
                b = mdl.Coefficients.Estimate; b = b(2:end);
            end
            if ir <= 3
                fit.Lambda1SE{1,ir}(subi,con) = fitinfo.Lambda1SE;
                fit.LambdaMinMSE{1,ir}(subi,con) = fitinfo.LambdaMinMSE;
                fit.BetaMinMSE{1,ir}(1:nfeatures,subi,con) = b(:,fitinfo.IndexMinMSE);
                fit.Beta1SE{1,ir}(1:nfeatures,subi,con) = b(:,fitinfo.Index1SE);
            else
                fit.BetaMinMSE{1,ir}(1:nfeatures,subi,con) = b;                
            end
        end
    end
    Nfeatures(con) = nfeatures;
end

varname = {'Raw data','Sample mean','Sample SD','Pop SD'} ;
varname = {'Raw data'} ;

save regression_result3_X1X30 Nfeatures fit varname



