function [samples] = ndRandn(meanmat,covmat,num)

if nargin == 2
    num = 1;
end

%identity: Y = MX --> I = MCxMt --> I = MVDVtMt --> MV(D^1/2)(D^1/2)VtMt -->
%M = (D^(-1/2)*V^T)
[V, D] = eig(covmat);
M = (sqrt(D)*V')';
X = randn(num,length(covmat));
samples = M * X';
samples = samples + meanmat ; 
end
