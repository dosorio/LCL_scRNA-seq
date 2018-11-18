function M=i_Magic_impute(data)

addpath('MAGIC_matlab\')
%sortby="expr_level";
%i_common_code;
%X=GM12878_expr(:,cellcycleGM12878=="G1");


% gene_names=cellstr(gl123);  MAGIC needs [cells x genes]
% data=X';

% library size normalization
libsize = sum(data,2);
data = bsxfun(@rdivide, data, libsize) * median(libsize);
% log transform -- usually one would log transform the data. Here we don't do it.
% data = log(data + 0.2);

[pc_imputed, U, ~] = run_magic(data, 'npca', 100, 'k', 15, 'a', 15, 'make_plot_opt_t', true);

% plot_genes = {'Cdh1', 'Vim', 'Fn1', 'Zeb1'};
% [M_imputed, genes_found] = project_genes(plot_genes, gene_names, pc_imputed, U);
M = pc_imputed * U';    % project

