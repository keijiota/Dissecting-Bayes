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

    for subi = 1:N
        endpoint = endpointorg(subi, :)'; 
        optaim = optAim(subi,:)';
        sample = ssamples(:,:,subi);
        PopSD  = orgnoise(subi,:); 
        
        if con == 1 || con == 2
           orgsample = sample - endpoint ;           
           X = [orgsample(:,1), orgsample(:,15), orgsample(:,30)];
        else
           orgsample = sample(:,1:5) - endpoint ;
           X = [orgsample(:,1), orgsample(:,3), orgsample(:,5)];
        end        
        % sorted raw data (not normalized)
        
        corrx = corr(X);
        CorrX(subi,1,con) = corrx(1,2);
        CorrX(subi,2,con) = corrx(1,3);
        CorrX(subi,3,con) = corrx(2,3);
        [R2 vif] = multicollinearity(X);
        VIF(subi,:,con) = vif;

        Y1 = optaim ; 
        Y2 = endpoint; 
        
        mdl1 = fitglm(X,Y1);
        mdl2 = fitglm(X,Y2);

        coefficient_opt(subi,:,con) = mdl1.Coefficients.Estimate(2:end); 
        coefficient_obs(subi,:,con) = mdl2.Coefficients.Estimate(2:end);
    end 
end
varname = {'Smin','S50%','Smax'} ; 
 
save regression_partial coefficient_opt coefficient_obs VIF varname









