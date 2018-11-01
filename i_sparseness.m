function s=i_sparseness(x)
n=numel(x);
s=(sqrt(n)-sum(abs(x))/norm(x))/(sqrt(n)-1);

% http://www.jmlr.org/papers/volume5/hoyer04a/hoyer04a.pdf 
% https://github.com/aludnam/MATLAB/tree/master/nmfpack
