load dataGSE65525
Xori=table2array(data(:,2:end));
scaleFactor = mean(sum(Xori))./sum(Xori);
X=Xori.*scaleFactor;

m=size(X,1);
xdata = mean(X,2);
cv=std(X,1,2)./xdata;
ydata = log10(cv);

xdata=xdata(~isnan(ydata));
ydata=ydata(~isnan(ydata));

xx=log10(xdata);

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
ysmooth = malowess(xx,ydata);
toc


[fitLoc]=fit(xx,ydata,'smoothingspline');

%%
figure;
scatter(xx,ydata);
hold on
scatter(xx,ysmooth)
% scatter(xx,fitLoc(xx))



%%
xSeq=min(xx):0.005:max(xx);
gapNum=zeros(length(xSeq),1);
for k=1:length(xSeq)
    gapNum(k)=sum(xx>= xSeq(k) - 0.05 & xx < xSeq(k) + 0.05);    
end

  cdx = find (gapNum > m*0.005);
  xSeq = 10.^xSeq;
  ySeq = fitLoc(log10(xSeq));

  %%
figure;
scatter(xx,ydata);
hold on
scatter(log10(xSeq),ySeq)
