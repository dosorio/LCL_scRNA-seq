% s4_hvgs;
% addpath('scGEApp\src');
datatag='../smpl0_GM12878_scRNA_seq_original';

load(datatag,'X0','core_idx0','genelist','s_phate0g1','cellcycletag0');
% X=X0(:,cellcycletag0=="G1");
% X=X0;
X=run_magic(X0,true);

i1=find(genelist=="REL");
i2=find(genelist=="RELA");
i3=find(genelist=="PRDM1");
i4=find(genelist=="AICDA");

i5=find(genelist=="PAX5");
i6=find(genelist=="BCL6");

i7=find(genelist=="IRF4");

i8=find(genelist=="BACH2");

x1=X(i1,:); x2=X(i2,:); x3=X(i3,:); x4=X(i4,:);
x5=X(i5,:); x6=X(i6,:); x7=X(i7,:); x8=X(i8,:);

% figure;
% subplot(2,2,1), stem(x1,'marker','none')
% subplot(2,2,2), stem(x2,'marker','none')
% subplot(2,2,3), stem(x3,'marker','none')
% subplot(2,2,4), stem(x4,'marker','none')
% 
% figure;
% subplot(2,2,1), stem(x5,'marker','none')
% subplot(2,2,2), stem(x6,'marker','none')
% subplot(2,2,3), stem(x7,'marker','none')
% subplot(2,2,4), stem(x8,'marker','none')

DataTable=table(x1',x2',x3',x4',x5',x6',x7',x8');
DataTable.Properties.VariableNames={'REL','RELA','PRDM1','AICDA','PAX5','BCL6','IRF4','BACH2'};

corrplot(DataTable)

%%
figure; scatter(x3,x1,'.'); xlabel('PRDM1 (Blimp)'); ylabel('REL (cRel)'); 
figure; scatter(x1,x4,'.'); xlabel('REL (cRel)'); ylabel('AICDA (Aicda)');
figure; scatter(x3,x4,'.'); xlabel('PRDM1 (Blimp)'); ylabel('AICDA (Aicda)');

%%
figure; gscatter(x3,x1,cellcycletag0); xlabel('PRDM1 (Blimp)'); ylabel('REL (cRel)'); 
figure; gscatter(x1,x4,cellcycletag0); xlabel('REL (cRel)'); ylabel('AICDA (Aicda)');
figure; gscatter(x3,x4,cellcycletag0); xlabel('PRDM1 (Blimp)'); ylabel('AICDA (Aicda)');

%%
figure; scatter3(x1,x3,x4,'.')
i=cellcycletag0=="G1";
figure; scatter3(x1(i),x3(i),x4(i),'.')
