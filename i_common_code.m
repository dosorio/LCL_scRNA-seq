load scExpr3GMs.mat

GM12878_expr=double(GM12878_expr);
GM18502_expr=double(GM18502_expr);
GMmix_expr=double(GMmix_expr);
n1=size(GM12878_expr,2);
n2=size(GM18502_expr,2);
n3=size(GMmix_expr,2);
[tagx1{1:n1}]=deal('eur');
[tagx2{1:n2}]=deal('afr');
[tagx3{1:n3}]=deal('mix');

switch sortby
    case "expr_level"
        S=sum(GM12878_expr,2)+sum(GM18502_expr,2)+sum(GMmix_expr,2);
        [expr_sorted,idx]=sort(S,'descend');
        GM12878_expr=GM12878_expr(idx,:);
        GM18502_expr=GM18502_expr(idx,:);
        GMmix_expr=GMmix_expr(idx,:);
        gl123=gl123(idx);
        gl123desc=gl123desc(idx);
    case "diff_expr"
        m_eur_gt_afr=(mean(GM12878_expr,2)-mean(GM18502_expr,2))>0;
        mdiff=abs(mean(GM12878_expr,2)-mean(GM18502_expr,2));
        mglob=mean([GM12878_expr GM18502_expr],2);
        [d_sorted,idx]=sort(mdiff./mglob,'descend');
        GM12878_expr=GM12878_expr(idx,:);
        GM18502_expr=GM18502_expr(idx,:);
        GMmix_expr=GMmix_expr(idx,:);
        gl123=gl123(idx);
        gl123desc=gl123desc(idx);
        mglob=mglob(idx);
        mdiff=mdiff(idx);
        m_eur_gt_afr=m_eur_gt_afr(idx);
       
    case "none"
        
    otherwise
        
end
%%
GMmix_afr=GMmix_expr(:,isafr==1);
GMmix_eur=GMmix_expr(:,isafr==-1);

%%
i1=(cellcycleGM12878=="G1");
i2=(cellcycleGM12878=="G2M");
i3=(cellcycleGM12878=="S");
j1=(cellcycleGM18502=="G1");
j2=(cellcycleGM18502=="G2M");
j3=(cellcycleGM18502=="S");

assert(size(GMmix_eur,2)==size(mixcellcycleGM12878,1))
assert(size(GMmix_afr,2)==size(mixcellcycleGM18502,1))

GMmix_afr_G1=GMmix_afr(:,mixcellcycleGM18502=="G1");
GMmix_afr_G2M=GMmix_afr(:,mixcellcycleGM18502=="G2M");
GMmix_afr_S=GMmix_afr(:,mixcellcycleGM18502=="S");

GMmix_eur_G1=GMmix_eur(:,mixcellcycleGM12878=="G1");
GMmix_eur_G2M=GMmix_eur(:,mixcellcycleGM12878=="G2M");
GMmix_eur_S=GMmix_eur(:,mixcellcycleGM12878=="S");



