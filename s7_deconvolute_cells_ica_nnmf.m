%% toy example
% S is matrix of expression levels of 5 genes x 2 cell types
S = [10.0000    1.0000
    9.0000    1.5000
    2.0000    8.0000
    1.0000    9.0000
    4.0000    5.0000];
% A is matrix of proportions of 2 cell types x 3 samples
A = [0.8000    0.9000    0.2000
    0.2000    0.1000    0.8000];

X=S*A;  % row = 5 genes, column = 3 cell types

[W,H]=nnmf(S*A,2)
    
%%
tic
[W,H]=nnmf(GMmix_expr,2);
toc

Hf=H./sum(H);
