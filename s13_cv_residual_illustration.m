sortby="expr_level";
i_common_code;


X=GM12878_expr(:,cellcycleGM12878=="G1");
% Y=GM18502_expr(:,cellcycleGM18502=="G1");
% Z=GM12878_expr(:,cellcycleGM12878=="S");
u=mean(X,2);
cv=std(X,0,2)./u;

%%
%close all
figure;
loglog(u,cv,'bo');
hold on
i=1e-4:0.4:1e2;
a=0.5;
j=10.^-(log10(i)*a);
plot(i,j,'-rv');
%j2=i.^-0.5;
%loglog(i,j2,'-o');
grid on
xlabel('Mean')
ylabel('CV')
% dt = datacursormode;
% dt.UpdateFcn = {@i_myupdatefcn,gl123};


%
% x = nbinrnd(20,.5,1000,1);
% params = nbinfit(x)                                   % unconstrained fit
% params = mle(x,'pdf',@(x,r)nbinpdf(x,r,.5),'start',23)     % constrain P=0.5
% params = mle(x,'pdf',@(x,r,p)nbinpdf(x,r,p),'start',[23 0.45])

% lgcv=log10(cv);
% 
% % lgcv ~ u
% % y=0.5*log10(b/x+a)
% 
% a=0.6028;
% b=1.084;
% 
% i=1e-4:1e3;
% j=10.^(0.5*log10(b./i+a));
% plot(i,j,'-r');

% a =       1.175;
% j=10.^(0.5*log10(b./i+a));
% plot(i,j,'-gx');


%

% [xData, yData] = prepareCurveData( u, cv );
% ft = fittype( '10.^(0.5*log10(b/x+a))', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = [0.6 1.84];
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );

lgcv=log10(cv);
[xData, yData] = prepareCurveData( u, lgcv );
ft = fittype( '0.5*log10(b/x+a)', 'independent', 'x', 'dependent', 'y' );
fo = fitoptions( 'Method', 'NonlinearLeastSquares' );
%opts.Display = 'Off';
%opts.StartPoint = [0.6 1.84];
% Fit model to data.
[fr] = fit( xData, yData, ft, fo );
ab=coeffvalues(fr);

i=1e-4:0.03:2.9e2;
j=10.^(0.5*log10(ab(2)./i+ab(1)));
plot(i,j,'-gs');

legend({'Genes','Poisson distribution',...
    sprintf('Nonlinear least squares\n(a=%.3f, b=%.3f)',...
    ab(1),ab(2))})

ylim([0.3 150])

% aabb=confint(fr);
% j=10.^(0.5*log10(aabb(1,1)./i+aabb(2,2)));
% plot(i,j,'--r');
% j=10.^(0.5*log10(aabb(2,1)./i+aabb(1,2)));
% plot(i,j,'--r');



