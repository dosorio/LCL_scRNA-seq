sortby="none";
i_common_code;

%% needs 3 mins each
% tic
% X=[GM12878_expr GM18502_expr GMmix_expr]';
% Y=tsne(X);
% toc
% Elapsed time is 2033.482995 seconds.

%%
% tic
% [~,s]=pca(X);
% figure;
% gscatter(s(:,1),s(:,2),tagx);
% toc
% % Elapsed time is 1985.682253 seconds.

% Y12878=tsne(GM12878_expr');
% Y18502=tsne(GM18502_expr');
% Ymix=tsne(GMmix_expr');
%%
load tsne_res % Y12878 Y18502 Ymix Y 

%%
figure;
subplot(2,2,1)
gscatter(Y(1:n1+n2,1),Y(1:n1+n2,2),[tagx1 tagx2]');
title('GM12878+GM18502')

subplot(2,2,2)
scatter(Y(1:n1+n2,1),Y(1:n1+n2,2),'w.');
hold on
scatter(Y(n1+n2+1:end,1),Y(n1+n2+1:end,2),'k.');
box on
title('GMmix')

subplot(2,2,3)
[~,i]=sort(cellcycleGM12878);
gscatter(Y12878(i,1),Y12878(i,2),cellcycleGM12878(i));
title('GM12878 (eur) only')
box on

subplot(2,2,4)
[~,i]=sort(cellcycleGM18502);
gscatter(Y18502(i,1),Y18502(i,2),cellcycleGM18502(i));
title('GM18502 (afr) only')
box on

% figure;
% scatter(Ymix(:,1),Ymix(:,2));
% title('GMmix')
% box on
