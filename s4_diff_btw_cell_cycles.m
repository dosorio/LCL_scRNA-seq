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

i1=(cellcycleGM12878=="G1");
i2=(cellcycleGM12878=="G2M");
i3=(cellcycleGM12878=="S");

j1=(cellcycleGM18502=="G1");
j2=(cellcycleGM18502=="G2M");
j3=(cellcycleGM18502=="S");

%%
close all
for k=101:110 % length(gl123)
    % if mean(GM12878_expr(k,:))<1, continue; end
    figure;
    subplot(2,2,1)
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
    subplot(2,2,2)
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
subplot(2,2,3)
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
end
