% load dataGSE65525
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

yy2 = smooth(log10(xdata),ydata,0.1,'rloess');

%%
figure;
scatter(log10(xdata),ydata);
hold on
scatter(log10(xdata),yy2)


xSeq=min(log10(xdata)):0.005:max(log10(xdata));
