load scExpr3GMs.mat
GM12878_expr=double(GM12878_expr);
GM18502_expr=double(GM18502_expr);
GMmix_expr=double(GMmix_expr);

%% needs 3 mins each
% Y12878=tsne(GM12878_expr');
% Y18502=tsne(GM18502_expr');
 Ymix=tsne(GMmix_expr');
%%
load tsne_res Y12878 Y18502 Ymix
%%
figure;
scatter(Y12878(:,1),Y12878(:,2),'.');
title('GM12878')
box on

figure;
scatter(Y18502(:,1),Y18502(:,2),'.');
title('GM18502')
box on

figure;
scatter(Ymix(:,1),Ymix(:,2),'.');
title('GMmix')
box on
