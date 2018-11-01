load dataGSE65525
Xori=table2array(data(:,2:end));
scaleFactor = mean(sum(Xori))./sum(Xori);
X=Xori.*scaleFactor;

m=size(X,1);
xdata = mean(X,2);
cv=std(X,1,2)./xdata;
lgydata = log10(cv);

xdata=xdata(~isnan(lgydata));

lgydata=lgydata(~isnan(lgydata));
lgxdata=log10(xdata);

%%
%  xdata <- xdata[is.na(ydata) != "TRUE"]
%  ydata <- ydata[is.na(ydata) != "TRUE"]
% fitLoc <- locfit.robust(ydata ~ lp(log10(xdata), nn = .2))

% takes 116 seconds
% tic
% ysmooth = smooth(xx,ydata,0.1,'rloess');
% toc


% takes 3 seconds
tic
ysmooth = malowess(lgxdata,lgydata);
toc

figure;
scatter(lgxdata,lgydata);
hold on
[~,idx]=sort(lgxdata);
plot(sort(lgxdata),ysmooth(idx),'-r');
% scatter(xx,ysmooth)
% scatter(xx,fitLoc(xx))

%%

[xData, yData] = prepareCurveData( 10.^lgxdata, lgydata );
ft = fittype( '0.5*log10(b/x+a)', 'independent', 'x', 'dependent', 'y' );
fo = fitoptions( 'Method', 'NonlinearLeastSquares' );
%opts.Display = 'Off';
%opts.StartPoint = [0.6 1.84];
% Fit model to data.
[fr] = fit( xData, yData, ft, fo );

ab=coeffvalues(fr);
%i=1e-4:0.03:2.9e2;
i1=-4:0.05:2;
i2=10.^i1;
j=(0.5*log10(ab(2)./i2+ab(1)));
plot(i1,j,'-g');


return;

%%

[fitLoc]=fit(xx,ydata,'poly2');
%[fitLoc]=fit(xx,ydata,'smoothingspline');
xSeq=min(xx):0.005:max(xx);
gapNum=zeros(length(xSeq),1);
for k=1:length(xSeq)
    gapNum(k)=sum(xx>= xSeq(k) - 0.05 & xx < xSeq(k) + 0.05);    
end

cdx = find (gapNum > m*0.005);
xSeq = 10.^xSeq;
ySeq = fitLoc(log10(xSeq));
  
figure;
scatter(xx,ydata);
hold on
scatter(log10(xSeq),ySeq)

%%
yDiff=diff(ySeq);
ix=find(yDiff>0 & xx(2:end)'>0);
if isempty(ix), ix=length(ySeq) - 1; end
    
xSeq_all=10.^(min(xx):0.001:max(xx));

  xSeq = xSeq(cdx(1):ix(1) + 1);
  ySeq = ySeq(cdx(1):ix(1) + 1);

  b =1;
  a =0;
  
  fun = @(x)0.5*log10(a(1) / x +a(2))-y;
  x = lsqnonlin(fun,[b a]);
  
  