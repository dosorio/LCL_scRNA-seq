sortby="diff_expr";
i_common_code;

%%
idx=mglob>=prctile(mglob,50);
GM12878_expr=GM12878_expr(idx,:);
GM18502_expr=GM18502_expr(idx,:);
GMmix_expr=GMmix_expr(idx,:);
gl123=gl123(idx);
gl123desc=gl123desc(idx);
mglob=mglob(idx);
mdiff=mdiff(idx);
m_eur_gt_afr=m_eur_gt_afr(idx);
d_sorted=d_sorted(idx);

T=table(gl123,gl123desc,mdiff,mglob,d_sorted,m_eur_gt_afr);

%%
close all
id=1316;
% figure;
% cdfplot(double(GM12878_expr(id,:)));
% hold on
% cdfplot(double(GM18502_expr(id,:)));
% legend({'GM12878','GM18502'})
% title(sprintf('%s [%s]',gl123(id),gl123desc(id)));

figure;
subplot(3,1,1)
stem(double(GM12878_expr(id,:)),'marker','none');
xlim([0 length(GM12878_expr(id,:))])
title('GM12878 EUR')

subplot(3,1,2)
stem(double(GMmix_expr(id,:)),'marker','none');
xlim([0 length(GMmix_expr(id,:))])
title('Mixture')


subplot(3,1,3)
stem(double(GM18502_expr(id,:)),'marker','none');
xlim([0 length(GM18502_expr(id,:))])
title('GM18502 AFR')
xlabel(sprintf('%s [%s]',gl123(id),gl123desc(id)));

