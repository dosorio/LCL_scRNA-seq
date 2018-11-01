sortby="expr_level";
i_common_code;


X=GM12878_expr(:,cellcycleGM12878=="G1");
Y=GM18502_expr(:,cellcycleGM18502=="G1");
% Z=GM12878_expr(:,cellcycleGM12878=="S");
u=mean(X,2);
cv=std(X,0,2)./u;

%%


figure;


sp=zeros(size(X,1),1);
for k=1:size(X,1)
    sp(k)=i_sparseness(X(k,:));
end

% idx=~isnan(cv);
% gl123x=gl123(idx);
% gl123descx=gl123desc(idx);
% u=u(idx);
% cv=cv(idx);

%%
figure;
loglog(u,sp,'o');

%%
figure;
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

%%

close all
u=mean(X,2);
cv=std(X,0,2)./u;
v=var(X,0,2);
i=u>0&v>0;

x=log10(u(i));
y=log10(v(i));
plot(x,y,'o');


p1 = polyfit(x,y,1);
p2 = polyfit(x,y,2);
p3 = polyfit(x,y,3);
p4 = polyfit(x,y,4);
p5 = polyfit(x,y,5);


v1 = polyval(p1,x);
v2 = polyval(p2,x);
v3 = polyval(p3,x);
v4 = polyval(p4,x);
v5 = polyval(p5,x);

y1=y./v1;
y2=y./v2;
y3=y./v3;
y4=y./v4;
y5=y./v5;

[r1,p1]=corr(x,y1,'type','k')
[r2,p2]=corr(x,y2,'type','k')
[r3,p3]=corr(x,y3,'type','k')
[r4,p4]=corr(x,y4,'type','k')
[r5,p5]=corr(x,y5,'type','k')




