function [res_cv]=i_logulogcv(C,rangev,uselog)

if nargin<3, uselog=true; end
if nargin<2 || isempty(rangev), rangev=[-5 3]; end

u=mean(C,2);
cv=std(C,[],2)./u;

lgcv=log10(cv);
[xData, yData] = prepareCurveData( u, lgcv );
ft = fittype( '0.5*log10(b/x+a)', 'independent', 'x', 'dependent', 'y' );
fo = fitoptions( 'Method', 'NonlinearLeastSquares' );
[fr] = fit( xData, yData, ft, fo );
ab=coeffvalues(fr);
if nargout>0
    lgcv_expected=(0.5*log10(ab(2)./u+ab(1)));
    res_cv=lgcv-lgcv_expected;
else
    if ~uselog
        scatter(log10(u),log10(cv))
        hold on
        fplot(@(x) -0.5*x, rangev,'r-')
        fplot(@(x) 0.5*log10(ab(2)./10.^x+ab(1)), rangev,'g-')
        box on
        xlabel('log10(u)');
        ylabel('log10(CV)');
    else
        loglog(u,cv,'o');
        hold on
        i=10^rangev(1):0.04:10^rangev(2);
        j=10.^-(0.5*log10(i));
        plot(i,j,'-rv');
        j=10.^(0.5*log10(ab(2)./i+ab(1)));
        plot(i,j,'-gs');
        grid on
        xlabel('u');
        ylabel('CV');
    end
end
