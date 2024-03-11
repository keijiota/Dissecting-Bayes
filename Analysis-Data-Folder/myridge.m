function [Bridge, fitinfo] = myridge(y,X,k,standardize,figureon)

% my ridge regression with cross validation and boot strapping
% y = ntrials * 1 matrix
% X = ntrials * nfeatures matrix
% k = k-fold cross validation, standardize = 0; figureon = 1;

ngrid = 100;
ncoeffs = size(X,2); ndata= size(X,1); nvalidate = ndata/k;
Xc = [ones(ndata,1),X];

iszero = 0; p = -1;
while ~ iszero % to determine the max lambda
    p = p + 1;
    lambda = 10^p;
    if standardize == 0
        Br = ridge(y,X,lambda,standardize);  % Regression includes the constant term
    elseif standardize == 1
        Br = ridge(y,Xc,lambda,standardize); % Regression doen't ibclude the constant term, thus add the constant term
    end
    iszero = abs(Br(2:end)) < 0.001 ;
    iszero = sum(iszero) == ncoeffs; % check if all the coefficients are 0
end
% lambda = logspace(-4,p,ngrid); % 10^4 ~ 10^p
maxlambda = 10^p; minlambda = 10^(-4);
a = logspace(2,p,round(ngrid*0.1));
b = logspace(-4,-1,round(ngrid*0.1));
c = linspace(max(b),min(a),ngrid-round(ngrid*0.1)*2);
lambda = [b,c,a];

bootN = 500; cvmse_boot = nan(k,ngrid,bootN); Br = nan(ncoeffs+1,ngrid,bootN);
for ib = 1:bootN
    indx  = randi(ndata,1,ndata);
    % indx  = [1:ndata];
    Br_boot = zeros(ncoeffs+1,ngrid);
    for ik = 1:k % cross validation
        tmp = ik*nvalidate-nvalidate+1:ik*nvalidate ;
        validate_set_ik = indx(tmp);
        train_set_ik = setdiff(indx,validate_set_ik);
        
        % split the data
        X_train   = X(train_set_ik,:);   X_validate  = X(validate_set_ik,:);
        Xc_train  = Xc(train_set_ik,:);  Xc_validate = Xc(validate_set_ik,:);
        y_train   = y(train_set_ik,1);   y_validate  = y(validate_set_ik,1);
        
        if standardize == 0
            br = ridge(y_train,X_train,lambda,standardize);
        elseif standardize == 1
            br = ridge(y_train,Xc_train,lambda,standardize);
        end
        y_hut = Xc_validate * br;
        
        Br_boot = Br_boot + br;
        cvmse_boot(ik,:,ib) = mean((y_validate - y_hut).^2);
    end
    Br(:,:,ib) = Br_boot / k;
end
% average across bootstrap samples
cvmse = mean(cvmse_boot,3); Bridge = mean(Br,3);
% mse across cross validation sets
mse = mean(cvmse); se = std(cvmse) ./ sqrt(k); 

[minMSE, indexminmse] = min(mse); lambdaMin = lambda(indexminmse);
minplus1 = mse(indexminmse) + se(indexminmse) ;

[~, index1se] = min((mse(indexminmse:end)-minplus1).^2) ;
index1se = indexminmse + index1se - 1;

if isempty(index1se)
    lambdaSE = [];
else
    lambdaSE = lambda(index1se);
end
% minMSE
% lambdaMin
% lambdaSE

fitinfo.intercept = Bridge(1,:);
Bridge = Bridge(2:end,:);

fitinfo.IndexMinMSE = indexminmse;
fitinfo.Index1SE = index1se;
fitinfo.MSE = mse;
fitinfo.minMSE = minMSE;
fitinfo.SE = se;
fitinfo.Lambda = lambda;
fitinfo.LambdaMinMSE = lambdaMin;
fitinfo.Lambda1SE = lambdaSE;

if figureon == 1
    figure;
    subplot(1,2,1); hold on
    errorbar(log(lambda), mse, se); xlim([log(minlambda) log(maxlambda)]);
    lineplot(log(lambdaMin),'v','k--'); lineplot(minplus1,'h','k--'); lineplot(log(lambdaSE),'v','b--');  ylabel('MSE');
    
    subplot(1,2,2);
    semilogx(lambda,Bridge); xlim([minlambda maxlambda]); ylabel('Beta');
end

