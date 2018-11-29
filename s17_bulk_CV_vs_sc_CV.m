% sortby="expr_level";
i_common_code;
i_get_genes_residual_cv_value;

clearvars -except gl123 res_cv
%%
load('\\wdxtba361-1\disk4t\1000GenomeRNAseq\expr\PeerRPKM_462_all_table.mat')
% X=X(:,~Ts.isyri2);
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

figure; scatter(res_cv,log10(bulk_cv))
corr(res_cv,log10(bulk_cv))
dt = datacursormode;
dt.UpdateFcn = {@myupdatefcn,gname};

%%
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