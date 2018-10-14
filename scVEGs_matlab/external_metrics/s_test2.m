load xxx GM12878G1valuesCV
vgene=GM12878G1valuesCV.geneavg;

%%
[cvv,gvv]=xlsread('cv_eur_afr.xlsx','Sheet1');
% targetglist=string(gvv(2:end,4));
targetglist=string(gvv(:,2));
targetgvalu=cvv(:,9);

i=targetgvalu>0 & targetglist~="";
targetglist=targetglist(i);
targetgvalu=targetgvalu(i);
%%
[isvgene]=ismember(targetglist,vgene);
y=log2(targetgvalu);
figure;
evqtlplot(y,double(isvgene))
[~,p]=kstest2(y(isvgene),y(~isvgene))

figure;
cdfplot(y(isvgene))
hold on
cdfplot(y(~isvgene))


y1=y(isvgene);
m1=median(y1);
y2=y(~isvgene);
m2=zeros(1000,1);
for k=1:1000
    y2x=y2(randperm(length(y2)));
    y2x=y2x(1:length(y1));
    m2(k)=median(y2x);
end
figure;
histogram(m2)
vline(m1,'r-')
sum(m2<=m1)/1000







