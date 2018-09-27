sortby="none";
i_common_code;

targetg=["ABCC1","ADRBK1","AQP5","CDKN2A","CYP11B2","CYP26B1","ENPP1","GLIS1","GNB3","HIF1A","IL4","IRS1","LBP","MEF2A","P2RY2","PAWR","PRKG2","TBXA2R","TDP1","THBS2","TMEM121","TNFRSF4","TP53"]';
targetg=string(textread('genelist_dna_repair.txt','%s'));
%%
close all
for k=1:length(gl123)
    if ~ismember(gl123(k),targetg), continue; end
    if median(GM12878_expr(k,:))<0.15, continue; end
    figure;
    subplot(2,2,1)
        hold on
        h1=cdfplot(GM12878_expr(k,:));
        h2=cdfplot(GM18502_expr(k,:));
        set(h1,'color','r','linewidth',2,'linestyle','-') 
        set(h2,'color','b','linewidth',2,'linestyle','-') 
        legend({'Eur','Afr'})
        [~,p]=kstest2(GM12878_expr(k,:)',GM18502_expr(k,:)');
        title(sprintf('%s (p=%.2e)',gl123(k),p));
        xlabel(sprintf('%s',gl123desc(k)));
        box on
    subplot(2,2,2)
        hold on
        h1=cdfplot(GM12878_expr(k,i1));
        h2=cdfplot(GM18502_expr(k,j1));
        set(h1,'color','r','linewidth',2,'linestyle','-') 
        set(h2,'color','b','linewidth',2,'linestyle','-') 
        legend({'Eur','Afr'})
        title(sprintf('G1',gl123(k)));
        box on
    subplot(2,2,3)
        hold on
        h1=cdfplot(GM12878_expr(k,i2));
        h2=cdfplot(GM18502_expr(k,j2));
        set(h1,'color','r','linewidth',2,'linestyle','-') 
        set(h2,'color','b','linewidth',2,'linestyle','-') 
        legend({'Eur','Afr'})
        title(sprintf('G2M',gl123(k)));
        box on
    subplot(2,2,4)
        hold on
        h1=cdfplot(GM12878_expr(k,i3));
        h2=cdfplot(GM18502_expr(k,j3));
        set(h1,'color','r','linewidth',2,'linestyle','-') 
        set(h2,'color','b','linewidth',2,'linestyle','-') 
        legend({'Eur','Afr'})
        title(sprintf('S',gl123(k)));
        box on
end