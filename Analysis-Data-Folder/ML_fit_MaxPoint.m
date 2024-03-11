
function [estimates, log_likelihood] = ML_fit_MaxPoint(Subdata, Smax)
% Subdata and S max = T trials * M conditions
% m1 = 30 dots 0 pena, m2 = 30 dots, -500 pena, m3 = 5 dots 0 pena, m4 = 30 dots, -500 pena
% Smax is the maximum point of sample

[T M] = size(Subdata); 
model = @SmaxStrategy;
start_point = rand(1, 3)*5;   % Random starting point for each parameter
estimates = fminsearchbnd3(model, start_point, [], []);

[log_likelihood] = SmaxStrategy(estimates);

    function [log_likelihood] = SmaxStrategy(params) 
        % y = (PB - bias) - Smax, get the maximum point to the subjective boundary
        % assuming that each subject has different bias for different penalty conditions
        
        bias_0 = params(1);
        bias_500 = params(2);         
        mysigma = params(3); 
        
        for con = 1:4
            if con == 1 || con == 3
                bias = bias_0; 
            else
                bias = bias_500; 
            end
            yhut(:,con) = (180 - bias) - Smax(:,con) ; 
        end
        
        diff = sum(sum((Subdata - yhut).^2));
        log_likelihood = -sum(-log(mysigma) - diff/(2*(mysigma.^2)));

    end
end


