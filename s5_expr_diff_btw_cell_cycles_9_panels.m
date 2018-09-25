sortby="expr_level";
sortby="diff_expr";
i_common_code;

%%
close all
for k=1201:1700 % length(gl123)
    k
    if mean(GM12878_expr(k,:))<0.3, continue; end
    figure;
    subplot(3,3,1)
        hold on
        h1=cdfplot(double(GM12878_expr(k,i1)));
        h2=cdfplot(double(GM12878_expr(k,i2)));
        h3=cdfplot(double(GM12878_expr(k,i3)));
        set(h1,'color','r','linewidth',2,'linestyle','-') 
        set(h2,'color','g','linewidth',2,'linestyle','-') 
        set(h3,'color','b','linewidth',2,'linestyle','-') 
        legend({'G1','G2M','S'})
        title(sprintf('Eur %s',gl123(k)));
        box on
    subplot(3,3,4)
        hold on
        h1=cdfplot(double(GM18502_expr(k,j1)));        
        h2=cdfplot(double(GM18502_expr(k,j2)));
        h3=cdfplot(double(GM18502_expr(k,j3)));
        set(h1,'color','r','linewidth',2,'linestyle','-') 
        set(h2,'color','g','linewidth',2,'linestyle','-') 
        set(h3,'color','b','linewidth',2,'linestyle','-')         
        legend({'G1','G2M','S'})
        title(sprintf('Afr %s',gl123(k)));
        box on
subplot(3,3,7)
        hold on
        h=cdfplot(double(GM12878_expr(k,i1)));
        set(h,'color','r','linewidth',1,'linestyle','-')        
        h=cdfplot(double(GM12878_expr(k,i2)));
        set(h,'color','g','linewidth',1,'linestyle','-')
        h=cdfplot(double(GM12878_expr(k,i3)));
        set(h,'color','b','linewidth',1,'linestyle','-')
        h=cdfplot(double(GM18502_expr(k,j1)));
        set(h,'color','r','linewidth',1,'linestyle','--')
        h=cdfplot(double(GM18502_expr(k,j2)));
        set(h,'color','g','linewidth',1,'linestyle','--')
        h=cdfplot(double(GM18502_expr(k,j3)));
        set(h,'color','b','linewidth',1,'linestyle','--')
        legend({'Eur-G1','Eur-G2M','Eur-S','Afr-G1','Afr-G2M','Afr-S'})
        title(sprintf('%s %s',gl123(k),gl123desc(k)));
        box on

        subplot(3,3,3)
        stem(GM12878_expr(k,:),'marker','none');
        xlim([0 length(GM12878_expr(k,:))])
        title('GM12878 EUR')

        subplot(3,3,6)
        stem(GMmix_expr(k,:),'marker','none');
        xlim([0 length(GMmix_expr(k,:))])
        title('Mixture')

        subplot(3,3,9)
        stem(GM18502_expr(k,:),'marker','none');
        xlim([0 length(GM18502_expr(k,:))])
        title('GM18502 AFR')
        xlabel(sprintf('%s [%s]',gl123(k),gl123desc(k)));       
        
        
        
    subplot(3,3,2)
        hold on
        h1=cdfplot(GM12878_expr(k,i1));
        h2=cdfplot(GM18502_expr(k,j1));
        set(h1,'color','r','linewidth',2,'linestyle','-') 
        set(h2,'color','b','linewidth',2,'linestyle','-') 
        legend({'Eur','Afr'})
        title(sprintf('G1',gl123(k)));
        box on
    subplot(3,3,5)
        hold on
        h1=cdfplot(GM12878_expr(k,i2));
        h2=cdfplot(GM18502_expr(k,j2));
        set(h1,'color','r','linewidth',2,'linestyle','-') 
        set(h2,'color','b','linewidth',2,'linestyle','-') 
        legend({'Eur','Afr'})
        title(sprintf('G2M',gl123(k)));
        box on
    subplot(3,3,8)
        hold on
        h1=cdfplot(GM12878_expr(k,i3));
        h2=cdfplot(GM18502_expr(k,j3));
        set(h1,'color','r','linewidth',2,'linestyle','-') 
        set(h2,'color','b','linewidth',2,'linestyle','-') 
        legend({'Eur','Afr'})
        title(sprintf('S',gl123(k)));
        box on
        
        
        
        
        
        
        
        
end
