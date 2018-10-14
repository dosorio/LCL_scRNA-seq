load xxx
vgene=GM12878G1valuesCV.geneavg;


%%
[y]=ismember(normexprlcl.gname,vgene);
sum(y)
x=table2array(normexprlcl(:,3:end));
x1=mean(x,2);
x2=std(x,0,2);
x3=x2./x1;

% figure;
% evqtlplot(x1,double(y))
% figure;
% evqtlplot(x2,double(y))
% figure;
% evqtlplot(x3,double(y))

%%
B=braingtexlassoh2(braingtexlassoh2.VarName4>0,:);
[isvgene]=ismember(B.XRRA1,vgene);
sum(isvgene)
y=B.VarName4;
figure;
evqtlplot(y,double(isvgene))
[~,p]=kstest2(y,double(isvgene))

figure;
cdfplot(y(isvgene))
hold on
cdfplot(y(~isvgene))

%%
H=H2Lassopredallgenes1mb500kb(H2Lassopredallgenes1mb500kb.Blood>0,:);

[isvgene]=ismember(H.c,vgene);
sum(isvgene)
y=H.Blood;
figure;
evqtlplot(y,double(isvgene))
[~,p]=kstest2(y,double(isvgene))


figure;
cdfplot(y(isvgene))
hold on
cdfplot(y(~isvgene))

%%
H=cveurg;
[isvgene]=ismember(H.A1BG,vgene);
sum(isvgene)
y=log(H.VarName2);
figure;
evqtlplot(y,double(isvgene))
[~,p]=kstest2(y(isvgene),y(~isvgene));

figure;
cdfplot(y(isvgene))
hold on
cdfplot(y(~isvgene))

%%
y1=y(isvgene);
m1=median(y1);
y2=y(~isvgene);

m2=zeros(10000,1);
for k=1:10000
    y2x=y2(randperm(length(y2)));
    y2x=y2x(1:length(y1));
    m2(k)=median(y2x);
end
figure;
histogram(m2)
vline(m1,'r-')








