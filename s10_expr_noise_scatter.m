sortby="expr_level";
i_common_code;

%%

X=GM12878_expr(:,cellcycleGM12878=="G1");
Y=GM18502_expr(:,cellcycleGM18502=="G1");
% Z=GM12878_expr(:,cellcycleGM12878=="S");

%%
v1=log10(std(X,0,2)./mean(X,2))-(-0.5*log10(mean(X,2)));
v2=log10(std(Y,0,2)./mean(Y,2))-(-0.5*log10(mean(Y,2)));
figure;
scatter(v1,v2,'.')
refline(1)
box on
xlabel('EUR')
ylabel('AFR')
dt = datacursormode;
dt.UpdateFcn = {@myupdatefcn,gl123};



%%

%{
figure;
u=mean(X,2);
cv=std(X,0,2)./u;

% idx=~isnan(cv);
% gl123=gl123(idx);
% gl123desc=gl123desc(idx);
% u=u(idx);
% cv=cv(idx);

loglog(u,cv,'.');
hold on
i=1e-4:1e2;
a=0.5;
j=10.^-(log10(i)*a);
plot(i,j,'-');

xlabel('Mean')
ylabel('CV')

% loglog(mean(Y,2),std(Y,[],2)./mean(Y,2),'.')
v=log10(cv)-(-0.5*log10(u));
%}



function txt = myupdatefcn(~,event_obj,g)
% Customizes text of data tips
% pos = event_obj.Position;
idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
txt = {g(idx)};
end