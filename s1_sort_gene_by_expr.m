load scExpr3GMs.mat
s=sum(GM12878_expr,2)+sum(GM18502_expr,2)+sum(GMmix_expr,2);
[s_sorted,idx]=sort(s,'descend');

%%
GM12878_expr=GM12878_expr(idx,:);
GM18502_expr=GM18502_expr(idx,:);
GMmix_expr=GMmix_expr(idx,:);
gl123=gl123(idx);
gl123desc=gl123desc(idx);
%%
id=503;
figure;
stem(sort(GM12878_expr(id,:),'descend'),'marker','none')
title(sprintf('%s [%s]',gl123(id),gl123desc(id)));
xlabel('Cell ID (sorted)')
ylabel('# of reads')



