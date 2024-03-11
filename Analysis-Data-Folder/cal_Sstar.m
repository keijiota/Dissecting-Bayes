function [S] = cal_Sstar(samplesigma,N,eg,prior_var,v_hut,Y)

% mearuement
sample_var = samplesigma^2;

% likelihood of population mean/sigma P(sample sigma | population sigma)
likelihood_var  = (N-1)/(2^((N-1)/2)*gamma((N-1)/2)*sample_var) .* ((v_hut*(N-1))/sample_var).^((N-1)/2-1) .* exp(-v_hut*(N-1)/(2*sample_var));  likelihood_var = likelihood_var/sum(likelihood_var);

% posterior of population mean/sigma P(population sigma | sample sigma)
posterior_var = prior_var .* likelihood_var; posterior_var = posterior_var/sum(posterior_var);

%---- calculation of expected gain -----
ExpectedExpectedGain = nan(length(Y), 2);

for iy = 1:length(Y)
    ExpectedExpectedGain(iy,:) = [nansum(eg(iy,:).*posterior_var), Y(iy)]; % expected expected gain taking a posterior into acount
end

[opteg tmp] = max(ExpectedExpectedGain(:,1));
S = ExpectedExpectedGain(tmp,2);

end