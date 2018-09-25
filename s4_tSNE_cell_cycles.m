sortby="none";
i_common_code;

%%

% X1=GM12878_expr(:,i1)';  % G1
% X2=GM12878_expr(:,i2)';  % G2M
% X3=GM12878_expr(:,i3)';  % S

X1=GM18502_expr(:,j1)';  % G1
X2=GM18502_expr(:,j2)';  % G2M
X3=GM18502_expr(:,j3)';  % S

Y1=tsne(X1);
Y2=tsne(X2);
Y3=tsne(X3);

%%
figure;
subplot(2,2,1)
scatter(Y1(:,1),Y1(:,2),'.');
title('G1')
box on
subplot(2,2,2)
scatter(Y2(:,1),Y2(:,2),'.');
title('G2M')
box on
subplot(2,2,3)
scatter(Y3(:,1),Y3(:,2),'.');
title('S')
box on


