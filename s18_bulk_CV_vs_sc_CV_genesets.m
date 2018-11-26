sortby="expr_level";
i_common_code;

% a=gl123; a(7568)="TT";
% i=find(extractBefore(a,3)~="IG");   % exclude all IG genes
% gl123=gl123(i);
% GM12878_expr=GM12878_expr(i,:);

X=GM12878_expr(:,cellcycleGM12878=="G1");
% Y=GM18502_expr(:,cellcycleGM18502=="G1");
% Z=GM12878_expr(:,cellcycleGM12878=="S");
u=mean(X,2);
cv=std(X,0,2)./u;

lgcv=log10(cv);
[xData, yData] = prepareCurveData( u, lgcv );
ft = fittype( '0.5*log10(b/x+a)', 'independent', 'x', 'dependent', 'y' );
fo = fitoptions( 'Method', 'NonlinearLeastSquares' );
[fr] = fit( xData, yData, ft, fo );
ab=coeffvalues(fr);
lgcv_expected=(0.5*log10(ab(2)./u+ab(1)));
res_cv=lgcv-lgcv_expected;
i=~isnan(res_cv);
gl123=gl123(i);
res_cv=res_cv(i);
X=X(i,:);
% u=u(i);
% figure; scatter(u(2:end),res_cv(2:end));

Trefg=readtable('autosomal_protein_coding_genes.txt');
i=ismember(gl123,Trefg.GENE);
gl123=gl123(i);
res_cv=res_cv(i);
X=X(i,:);

% clearvars -except gl123 res_cv
%%
load('\\wdxtba361-1\disk4t\1000GenomeRNAseq\expr\PeerRPKM_462_all_table.mat')
%X=X(:,~Ts.isyri2);
X=X(:,Ts.popid2=="CEU");
u=mean(X,2);
bulk_cv=std(X,0,2)./u;
i=bulk_cv>0;
bulk_cv=bulk_cv(i);
Tg=Tg(i,:);

[gname,i,j]=intersect(gl123,Tg.genename);
res_cv=res_cv(i);
bulk_cv=bulk_cv(j);
u=u(j);

% figure; scatter(res_cv,bulk_cv)
% corr(res_cv,bulk_cv)

clearvars -except res_cv bulk_cv gname
%%
% load('\\cvm-research-dr.cvm.tamu.edu\cailab\Cai-Terra3\DATA\DISK4T\ref_gene_sets\msigdb_v61_c7.mat');
load('\\cvm-research-dr.cvm.tamu.edu\cailab\Cai-Terra3\DATA\DISK4T\ref_gene_sets\msigdb_v62\msigdb_v62_c5.mat');
T=cell2table(cell(0,8),'VariableNames',{'SetID','R1','P1','R2','P2','R3','P3','SetName'});
c=1;
for k=1:length(GeneSet)
    [~,idx]=intersect(gname,string(GeneSet{k})');
    if numel(idx)<5, continue; end
    [r1,p1]=corr(res_cv(idx),log10(bulk_cv(idx)));
    [r2,p2]=corr(res_cv(idx),log10(bulk_cv(idx)),'type','s');
    [r3,p3]=corr(res_cv(idx),log10(bulk_cv(idx)),'type','k');
    if any([p1 p2 p3]<0.01) && any([abs(r1) abs(r2) abs(r3)]>0.3)
        T=[T;{k,r1,p1,r2,p2,r3,p3,GeneSetName{k}}];
%         T.id(c)=k;
%         T.r1(c)=r1; T.p1(c)=p1;
%         T.r2(c)=r2; T.p2(c)=p2;
%         T.r3(c)=r3; T.p3(c)=p3;  
%         T.gsetname{c}=GeneSetName{k};
        c=c+1;
    end
end

T2=T(T.P3<0.05/1000,:);
T2=sortrows(T2,6,'descend');


%%
close all
for k=1:9
    targetk=T2.SetID(k);
    [gnamex,idx]=intersect(gname,string(GeneSet{targetk}));
    figure;
    scatter(res_cv(idx),log10(bulk_cv(idx)),'ro');
    refline
    text(0.01+res_cv(idx),log10(bulk_cv(idx)),gnamex,'rotation',0);   

    title(strrep(GeneSetName{targetk},'_','\_'));
    [r1,p1]=corr(res_cv(idx),log10(bulk_cv(idx)));
    [r2,p2]=corr(res_cv(idx),log10(bulk_cv(idx)),'type','s');
    [r3,p3]=corr(res_cv(idx),log10(bulk_cv(idx)),'type','k');
    xlabel('Residual CV')
    ylabel('Population CV')
    % dt = datacursormode;
    % dt.UpdateFcn = {@myupdatefcn,gnamex};
    box on
end



function txt = myupdatefcn(~,event_obj,g)
% Customizes text of data tips
pos = event_obj.Position;
idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
num2str(pos(1))
num2str(pos(2))
txt = {char(g(idx))};
% txt={num2str(pos(2))}
end