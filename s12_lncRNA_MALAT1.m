% We propose that non-linear sequence homology—in which the relative abundance of a set of protein-binding motifs is conserved, but the sequential relationships between them are not—is prevalent in lncRNAs. To quantify non-linear homology, we introduce SEEKR,

% TESTING: lncRNA with non-linear homology should co-expressed


% https://www.nature.com/articles/s41588-018-0207-8
% Unexpected relationships also emerged in the hierarchical clusters 
% of Fig. 2. Most notably, the lncRNAs NEAT1 and MALAT1 showed greater 
% than average similarity to XIST in both human and mouse.

i1=find(gl123=="MALAT1");
i2=find(gl123=="NEAT1");
i3=find(gl123=="XIST");

i1=i3;

figure;
subplot(2,2,1)
scatter(GM12878_expr(i1,:),GM12878_expr(i2,:))
subplot(2,2,2)
scatter(GM18502_expr(i1,:),GM18502_expr(i2,:))

figure;
subplot(2,1,1)
plot(GM12878_expr(i1,:))
subplot(2,1,2)
plot(GM12878_expr(i2,:))

figure;
subplot(2,1,1)
plot(GM18502_expr(i1,:))
subplot(2,1,2)
plot(GM18502_expr(i2,:))

% https://www.nature.com/articles/s41588-018-0252-3
% Notably, this Malat1-deletion model exhibited substantial upregulation
% of Malat1’s adjacent genes, including Neat1, Frmd8, Tigd3, Ehbp1l1, 
% Ltbp3,

