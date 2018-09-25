%sortby="expr_level";
sortby="diff_expr";
%sortby="none";
i_common_code;

%%
close all
% c=0;
% for k=1:length(gl123)
%     k
%     if median(GM12878_expr(k,:))<0.25, continue; end
%     c=c+1;
% end
% c

%%
for k=1:length(gl123)
    k
    if median(GM12878_expr(k,:))<0.25, continue; end
    
    fh=figure('visible','off');
    subplot(3,2,1)
        hold on
        h1=cdfplot(double(GM12878_expr(k,i1)));
        h2=cdfplot(double(GM12878_expr(k,i2)));
        h3=cdfplot(double(GM12878_expr(k,i3)));
        set(h1,'color',[102 194 164]./255,'linewidth',2,'linestyle','-') 
        set(h2,'color',[35,139,69]./255,'linewidth',2,'linestyle','-') 
        set(h3,'color',[0,68,27]./255,'linewidth',2,'linestyle','-') 
        legend({'G1','G2M','S'})
        title(sprintf('GM12878 EUR %s',gl123(k)));
        box on
        xlabel('')
    subplot(3,2,5)
        hold on
        h1=cdfplot(double(GM18502_expr(k,j1)));        
        h2=cdfplot(double(GM18502_expr(k,j2)));
        h3=cdfplot(double(GM18502_expr(k,j3)));
        set(h1,'color',[252,141,89]./255,'linewidth',2,'linestyle','-') 
        set(h2,'color',[215,48,31]./255,'linewidth',2,'linestyle','-') 
        set(h3,'color',[127,0,0]./255,'linewidth',2,'linestyle','-') 
        legend({'G1','G2M','S'})
        title(sprintf('GM18502 AFR %s',gl123(k)));
        box on
        xlabel('')
subplot(3,2,3)
        hold on
        h1=cdfplot(double(GM12878_expr(k,i1)));        
        h2=cdfplot(double(GM12878_expr(k,i2)));        
        h3=cdfplot(double(GM12878_expr(k,i3)));
        set(h1,'color',[102 194 164]./255,'linewidth',2,'linestyle','-') 
        set(h2,'color',[35,139,69]./255,'linewidth',2,'linestyle','-') 
        set(h3,'color',[0,68,27]./255,'linewidth',2,'linestyle','-') 
        
        h1=cdfplot(double(GM18502_expr(k,j1)));        
        h2=cdfplot(double(GM18502_expr(k,j2)));        
        h3=cdfplot(double(GM18502_expr(k,j3)));
        set(h1,'color',[252,141,89]./255,'linewidth',2,'linestyle','-') 
        set(h2,'color',[215,48,31]./255,'linewidth',2,'linestyle','-') 
        set(h3,'color',[127,0,0]./255,'linewidth',2,'linestyle','-') 
        
        % legend({'Eur-G1','Eur-G2M','Eur-S','Afr-G1','Afr-G2M','Afr-S'})
        % title(sprintf('%s %s',gl123(k),gl123desc(k)));
        % title(sprintf('%s',gl123(k)));
        title(sprintf('GM12878 & GM18502 %s',gl123(k)));
        box on
        xlabel('')

subplot(3,2,2)
        stem(GM12878_expr(k,:),'marker','none','color',[35,139,69]./255);
        xlim([0 length(GM12878_expr(k,:))])
        title(sprintf('GM12878 EUR %s',gl123(k)));

        subplot(3,2,4)
        stem(GMmix_expr(k,:),'marker','none','color',[115,115,115]./255);
        xlim([0 length(GMmix_expr(k,:))])
        title(sprintf('GM12878+GM18502 MIX %s',gl123(k)));

        subplot(3,2,6)
        stem(GM18502_expr(k,:),'marker','none','color',[215,48,31]./255);
        xlim([0 length(GM18502_expr(k,:))])
        title(sprintf('GM18502 AFR %s',gl123(k)));
        %xlabel(sprintf('%s [%s]',gl123(k),gl123desc(k)));
        xlabel(sprintf('%s [%s]',gl123(k),gl123desc(k)));
        if m_eur_gt_afr(k)
            print(fh,sprintf('img/eur/%s.png',gl123(k)),'-dpng');
        else
            print(fh,sprintf('img/afr/%s.png',gl123(k)),'-dpng');
        end
        %saveas(fh,sprintf('img/%s.png',gl123(k)));
        close(fh);        
end
