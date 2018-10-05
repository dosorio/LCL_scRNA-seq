sortby="expr_level";
i_common_code;

%%
% close all
% A=GM12878_expr;
A=GM12878_expr; a=cellcycleGM12878;
% A=GM18502_expr; a=cellcycleGM18502;
[x,i]=sort(a);
figure;
imagesc(log10(1+A(1:8000,i)))
vline(find(x=="G2M", 1 ),'y-')
vline(find(x=="S", 1 ),'g-')
title('G1 |  G2/M | S')
colorbar;
xlabel('Cell')
ylabel('Gene')

