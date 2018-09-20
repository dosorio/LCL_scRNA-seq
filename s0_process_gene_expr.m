gl1=string(textread('GM12878_gene.txt','%s'));
gl2=string(textread('GM18502_gene.txt','%s'));
gl3=string(textread('GMmix_gene.txt','%s'));

gl123=intersect(gl1,intersect(gl2,gl3));

[~,~,j1]=intersect(gl123,gl1);
[~,~,j2]=intersect(gl123,gl2);
[~,~,j3]=intersect(gl123,gl3);

%%
load GM12878_expr.txt
GM12878_expr=GM12878_expr(j1,:);
load GM18502_expr.txt
GM18502_expr=GM18502_expr(j2,:);
load GMmix_expr.txt
GMmix_expr=GMmix_expr(j3,:);

GM12878_expr=uint16(GM12878_expr);
GM18502_expr=uint16(GM18502_expr);
GMmix_expr=uint16(GMmix_expr);
%%
save scExpr3GMs GM12878_expr GM18502_expr GMmix_expr gl123




