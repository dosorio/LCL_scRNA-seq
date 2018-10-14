sortby="expr_level";
cd ..
i_common_code;
cd vs_pooled_rnaseq\

%%

ethid="EUR";
switch ethid
    case "EUR"
X=GM12878_expr; x=cellcycleGM12878; x1=i1; x2=i2; x3=i3;
    case "AFR"
X=GM18502_expr; x=cellcycleGM18502; x1=j1; x2=j2; x3=j3;
end

% X=X(:,i3);
scExprsum=sum(X,2);



%%
[cvv,gvv]=xlsread('cv_eur_afr.xlsx','Sheet1');
% targetglist=string(gvv(2:end,4));
targetglist=string(gvv(:,2));
targetgvalu=cvv(:,1);

%%
load ../scVEGs_matlab/external_metrics/xxx.mat isvgene

[glist,i,j]=intersect(gl123,targetglist);
T=table(glist, scExprsum(i), targetgvalu(j), isvgene(i));


%%
figure;
i=T.Var2>0&T.Var3>0;
T=T(i,:);
x=log(T.Var2);
y=log(T.Var3);

scatter(x,y,'.')
corr(x,y)

pause
% figure;
% scatter(x(T.Var4),y(T.Var4),'.')
% corr(x(T.Var4),y(T.Var4))

%%
% T.glist
% T.Var2
load gd660
[glist,i,j]=intersect(gl123,glist660);
T=table(glist, scExprsum(i), mn(j), isvgene(i));

i=T.Var2>0&T.Var3>0;
T=T(i,:);
x=log(T.Var2);
y=log(T.Var3);

scatter(x,y,'.')
corr(x,y)
