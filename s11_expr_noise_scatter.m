sortby="expr_level";
i_common_code;


X=GM12878_expr(:,cellcycleGM12878=="G1");
Y=GM18502_expr(:,cellcycleGM18502=="G1");
% Z=GM12878_expr(:,cellcycleGM12878=="S");

%%


figure;
u=mean(X,2);
cv=std(X,0,2)./u;

% idx=~isnan(cv);
% gl123x=gl123(idx);
% gl123descx=gl123desc(idx);
% u=u(idx);
% cv=cv(idx);

loglog(u,cv,'o');
hold on
i=1e-4:1e2;
a=0.5;
j=10.^-(log10(i)*a);
plot(i,j,'-rs');
%j2=i.^-0.5;
%loglog(i,j2,'-o');
grid on
xlabel('Mean')
ylabel('CV')
legend({'genes','Poisson'})

% loglog(mean(Y,2),std(Y,[],2)./mean(Y,2),'.')
% v=log10(cv)-(-0.5*log10(u));

%{
For single-cell gene expression data, in the
ideal condition all genes should obey CV = u^1/2 [11],
following a Poisson distribution as depicted by a black
diagonal line in log(?) vs log(CV) plot shown in Fig. 2.
BMC Genomics 2016, 17(Suppl 7):508

%}

%%
%{
v1=log10(std(X,0,2)./mean(X,2))-(-0.5*log10(mean(X,2)));
v2=log10(std(Y,0,2)./mean(Y,2))-(-0.5*log10(mean(Y,2)));

figure;
scatter(v1,v2,'.')
refline(1)
box on
xlabel('EUR residual cv')
ylabel('AFR residual cv')
dt = datacursormode;
dt.UpdateFcn = {@i_myupdatefcn,gl123};
%}

