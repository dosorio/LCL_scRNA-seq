% sortby="expr_level";
% i_common_code;

% a=gl123; a(7568)="TT";
% i=find(extractBefore(a,3)~="IG");   % exclude all IG genes
% gl123=gl123(i);
% GM12878_expr=GM12878_expr(i,:);

X=GM12878_expr(:,cellcycleGM12878=="G1");
% Y=GM18502_expr(:,cellcycleGM18502=="G1");
% Z=GM12878_expr(:,cellcycleGM12878=="S");
u=mean(X,2);
cv=std(X,0,2)./u;

lgcv=log10(cv);
[xData, yData] = prepareCurveData( u, lgcv );
ft = fittype( '0.5*log10(b/x+a)', 'independent', 'x', 'dependent', 'y' );
fo = fitoptions( 'Method', 'NonlinearLeastSquares' );
[fr] = fit( xData, yData, ft, fo );
ab=coeffvalues(fr);
lgcv_expected=(0.5*log10(ab(2)./u+ab(1)));
res_cv=lgcv-lgcv_expected;
i=~isnan(res_cv);
gl123=gl123(i);
res_cv=res_cv(i);
% u=u(i);
% figure; scatter(u(2:end),res_cv(2:end));

Trefg=readtable('autosomal_protein_coding_genes.txt');
i=ismember(gl123,Trefg.GENE);
gl123=gl123(i);
res_cv=res_cv(i);

% clearvars -except gl123 res_cv