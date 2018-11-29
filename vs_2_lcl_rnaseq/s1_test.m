cd ..
sortby="expr_level";
i_common_code;
cd vs_2_lcl_rnaseq\
%%

C1=GM12878_expr; %(:,cellcycleGM12878=="G1");
C2=GM18502_expr; %(:,cellcycleGM18502=="G1");
m1=mean(C1,2)*100;
m2=mean(C2,2)*100;
rmdv=(m1-m2)./mean([m1 m2],2);
Tsc=table(gl123,m1,m2,rmdv);

%%
Tsp=readtable('RNAseq.tsv.txt','ReadVariableNames',true);
Tsp.RMD=(Tsp.GM12878-Tsp.GM18502)./mean([Tsp.GM12878 Tsp.GM18502],2);

[glist,i,j]=intersect(gl123,Tsp.Gene,'stable');
Tsp=Tsp(j,:);
Tsc=Tsc(i,:);

%%
corrmyown(Tsc.m1,Tsp.GM12878);
corrmyown(Tsc.m2,Tsp.GM18502);
corrmyown(Tsc.m2,Tsp.GM12878);
corrmyown(Tsc.m1,Tsp.GM18502);


figure;
subplot(2,2,1); loglog(Tsc.m1,Tsp.GM12878,'o')
subplot(2,2,2); loglog(Tsc.m2,Tsp.GM18502,'o')
subplot(2,2,3); loglog(Tsc.m1,Tsp.GM18502,'o')
subplot(2,2,4); loglog(Tsc.m2,Tsp.GM12878,'o')
