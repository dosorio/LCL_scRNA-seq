sortby="expr_level";
i_common_code;

%%
%    case "EUR"
X=GM12878_expr; x=cellcycleGM12878; x1=i1; x2=i2; x3=i3;
A1=X(:,x1); [~,i]=sort(sum(A1,1),'descend'); A1=A1(:,i);
A2=X(:,x2); [~,i]=sort(sum(A2,1),'descend'); A2=A2(:,i);
A3=X(:,x3); [~,i]=sort(sum(A3,1),'descend'); A3=A3(:,i);
%    case "AFR"
X=GM18502_expr; x=cellcycleGM18502; x1=j1; x2=j2; x3=j3;
B1=X(:,x1); [~,i]=sort(sum(B1,1),'descend'); B1=B1(:,i);
B2=X(:,x2); [~,i]=sort(sum(B2,1),'descend'); B2=B2(:,i);
B3=X(:,x3); [~,i]=sort(sum(B3,1),'descend'); B3=B3(:,i);

%%
G1=GMmix_eur_G1; S=GMmix_eur_S; G2M=GMmix_eur_G2M;
sG1=sum(G1); sS=sum(S); sG2M=sum(G2M);
a=[sG1, sS, sG2M]';
b=[ones(size(sG1)), 2*ones(size(sS)), 3*ones(size(sG2M))]';

G1=GMmix_afr_G1; S=GMmix_afr_S; G2M=GMmix_afr_G2M;
sG1=sum(G1); sS=sum(S); sG2M=sum(G2M);
c=[sG1, sS, sG2M]';
d=[4*ones(size(sG1)), 5*ones(size(sS)), 6*ones(size(sG2M))]';
figure; boxplot(log10([a;c]),[b;d])
figure; cdfplot(a(b==1)); hold on; cdfplot(c(d==4))
[~,p]=kstest2(a(b==2),c(d==5))
figure; cdfplot(a(b==2)); hold on; cdfplot(c(d==5))
[~,p]=kstest2(a(b==2),c(d==5))
figure; cdfplot(a(b==3)); hold on; cdfplot(c(d==6))
[~,p]=kstest2(a(b==3),c(d==6))

%%
G1=GMmix_eur_G1; S=GMmix_eur_S; G2M=GMmix_eur_G2M;
sG1=sum(G1); sS=sum(S); sG2M=sum(G2M);
a=[sG1, sS, sG2M]';
b=[ones(size(sG1)), 2*ones(size(sS)), 3*ones(size(sG2M))]';

X=GM12878_expr; x=cellcycleGM12878; x1=i1; x2=i2; x3=i3;
A1=X(:,x1); [~,i]=sort(sum(A1,1),'descend'); A1=A1(:,i);
A2=X(:,x2); [~,i]=sort(sum(A2,1),'descend'); A2=A2(:,i);
A3=X(:,x3); [~,i]=sort(sum(A3,1),'descend'); A3=A3(:,i);

G1=A1; S=A3; G2M=A2;
sG1=sum(G1); sS=sum(S); sG2M=sum(G2M);
c=[sG1, sS, sG2M]';
d=[4*ones(size(sG1)), 5*ones(size(sS)), 6*ones(size(sG2M))]';
close all
figure; boxplot(log10([a;c]),[b;d])
figure; cdfplot(a(b==1)); hold on; cdfplot(c(d==4))
[~,p]=kstest2(a(b==1),c(d==4))
figure; cdfplot(a(b==2)); hold on; cdfplot(c(d==5))
[~,p]=kstest2(a(b==2),c(d==5))
figure; cdfplot(a(b==3)); hold on; cdfplot(c(d==6))
[~,p]=kstest2(a(b==3),c(d==6))

%%

G1=GMmix_afr_G1; S=GMmix_afr_S; G2M=GMmix_afr_G2M;
sG1=sum(G1); sS=sum(S); sG2M=sum(G2M);
a=[sG1, sS, sG2M]';
b=[ones(size(sG1)), 2*ones(size(sS)), 3*ones(size(sG2M))]';

X=GM18502_expr; x=cellcycleGM18502; x1=j1; x2=j2; x3=j3;
B1=X(:,x1); [~,i]=sort(sum(B1,1),'descend'); B1=B1(:,i);
B2=X(:,x2); [~,i]=sort(sum(B2,1),'descend'); B2=B2(:,i);
B3=X(:,x3); [~,i]=sort(sum(B3,1),'descend'); B3=B3(:,i);

G1=B1; S=B3; G2M=B2;
sG1=sum(G1); sS=sum(S); sG2M=sum(G2M);
c=[sG1, sS, sG2M]';
d=[4*ones(size(sG1)), 5*ones(size(sS)), 6*ones(size(sG2M))]';
close all
figure; boxplot(log10([a;c]),[b;d])
figure; cdfplot(a(b==1)); hold on; cdfplot(c(d==4))
[~,p]=kstest2(a(b==1),c(d==4))
figure; cdfplot(a(b==2)); hold on; cdfplot(c(d==5))
[~,p]=kstest2(a(b==2),c(d==5))
figure; cdfplot(a(b==3)); hold on; cdfplot(c(d==6))
[~,p]=kstest2(a(b==3),c(d==6))

