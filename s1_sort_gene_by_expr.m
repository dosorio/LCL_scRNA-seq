load scExpr3GMs.mat

%%
S=sum(GM12878_expr,2)+sum(GM18502_expr,2)+sum(GMmix_expr,2);
[expr_sorted,idx]=sort(S,'descend');
GM12878_expr=GM12878_expr(idx,:);
GM18502_expr=GM18502_expr(idx,:);
GMmix_expr=GMmix_expr(idx,:);
gl123=gl123(idx);
gl123desc=gl123desc(idx);

%%
median_eur=mean(GM12878_expr,2);
median_afr=mean(GM18502_expr,2);
median_mix=mean(GMmix_expr,2);

n=length(gl123);
pvalue=ones(n,1);
tic
for k=1:n
    [~,p]=kstest2(GM12878_expr(k,:),GM18502_expr(k,:));
    pvalue(k)=p;
end
toc
%[FDR,qvalue]=mafdr(pvalue);
%tag05=FDR<0.05;
bonferroni=pvalue<0.05/n;
T=table(gl123,gl123desc,median_eur,...
    median_mix,median_afr,pvalue,bonferroni);
%%

id=781;
figure;
subplot(3,2,[1 3 5])
% stem(sort(GM12878_expr(id,:),'descend'),'marker','none')
plot(sort(GM12878_expr(id,1:5530),'descend'),'r-','linewidth',2)
hold on
plot(sort(GMmix_expr(id,1:5530),'descend'),'g-','linewidth',2)
plot(sort(GM18502_expr(id,1:5530),'descend'),'b-','linewidth',2)
legend({'GM12878 (EUR)','Mixture','GM18502 (AFR)'})
xlim([-30 5540])
title(sprintf('%s [%s]',gl123(id),gl123desc(id)));
[~,p]=kstest2(GM12878_expr(id,:),GM18502_expr(id,:));
xlabel(sprintf('Cell ID (sorted), P=%.2e',...
      T.pvalue(id)))
ylabel('# of reads')


%figure;
subplot(3,2,2)

stem(double(GM12878_expr(id,:)),'marker','none','color','k');
xlim([0 length(GM12878_expr(id,:))])
title('GM12878 EUR')

subplot(3,2,4)
stem(double(GMmix_expr(id,:)),'marker','none','color','k');
xlim([0 length(GMmix_expr(id,:))])
title('Mixture')


subplot(3,2,6)
%subplot(3,1,3)
stem(double(GM18502_expr(id,:)),'marker','none','color','k');
xlim([0 length(GM18502_expr(id,:))])
title('GM18502 AFR')
xlabel(sprintf('%s [%s]',gl123(id),gl123desc(id)));


