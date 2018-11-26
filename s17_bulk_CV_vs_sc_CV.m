sortby="expr_level";
i_common_code;

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

clearvars -except gl123 res_cv
%%
load('\\wdxtba361-1\disk4t\1000GenomeRNAseq\expr\PeerRPKM_462_all_table.mat')
% X=X(:,~Ts.isyri2);
X=X(:,Ts.popid2=="CEU");
u=mean(X,2);
bulk_cv=std(X,0,2)./u;
i=bulk_cv>0;
bulk_cv=bulk_cv(i);
Tg=Tg(i,:);

[gname,i,j]=intersect(gl123,Tg.genename);
res_cv=res_cv(i);
bulk_cv=bulk_cv(j);
u=u(j);

% figure; scatter(res_cv,bulk_cv)
% corr(res_cv,bulk_cv)

figure; scatter(res_cv,log10(bulk_cv))
corr(res_cv,log10(bulk_cv))
dt = datacursormode;
dt.UpdateFcn = {@myupdatefcn,gname};

%%
function txt = myupdatefcn(~,event_obj,g)
% Customizes text of data tips
pos = event_obj.Position;
idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
num2str(pos(1))
num2str(pos(2))
txt = {char(g(idx))};
% txt={num2str(pos(2))}
end