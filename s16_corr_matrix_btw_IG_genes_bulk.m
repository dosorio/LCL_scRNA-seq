load('\\wdxtba361-1\disk4t\1000GenomeRNAseq\expr\PeerRPKM_462_all_table.mat')
% Y=X(:,Ts.isyri2);
Y=X(:,Ts.popid2=="CEU");
Z=Y(strlength(Tg.genename)>3,:);
Tg=Tg(strlength(Tg.genename)>3,:);
A=Z(extractBefore(Tg.genename,3)=="IG",:);
Tg=Tg(extractBefore(Tg.genename,3)=="IG",:);

g=Tg.genename;


load targetIGgenes target_ig_genes
i=ismember(g,target_ig_genes);
g=g(i);
A=A(i,:);
%%   
B=corr(A');
C=B>0.5;
C=C-diag(diag(C));
% C=logical(C);
% B(~C)=0;

i=sum(C)==0;
C(i,:)=[]; C(:,i)=[];
g2=g; g2(i)=[];

%%
% [x,y,v]=find(triu(B));
% for k=1:length(x)
%     
% end
bg1 = biograph(C>0,g2);
figure;
view(bg1);

%%
figure;
imagesc(B)
dt = datacursormode;
dt.UpdateFcn = {@myupdatefcn,g};
colorbar

% EUR modules:
% IGHA1 – IGLC2 – IGLV2-23
% IGHG1 – IGKC – IGHG3 – IGHG4 – IGKV4-1

% AFR modules:
% IGHG1 -  IGLC3 – IGLC2
% IGKV1-12 – IGKC – IGHM – IGKV3-15 – IGHV3-74


%{
%% Plot component scores.
%Create a plot of the first two columns of score.

% PCA needs [samples x features] here is [genes x cells]

[coeff,score,latent,tsquared,explained] = pca(A); % ,'VariableWeights','variance');
figure
plot(score(:,1),score(:,2),'+')
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')
xlabel(sprintf('PCA axis 1 (%.2f%%)',explained(1)))
ylabel(sprintf('PCA axis 2 (%.2f%%)',explained(2)))
for k=1:length(g)
    text(score(k,1),score(k,2),g(k));
end

% tSNE needs [samples x features] here is [genes x cells]
figure;
score=tsne(A);
scatter(score(:,1),score(:,2))
for k=1:length(g)
    text(score(k,1),score(k,2),g(k));
end
%%
%}


function txt = myupdatefcn(~,event_obj,g)
% Customizes text of data tips
pos = event_obj.Position;
% idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
num2str(pos(2))
txt = {[char(g(pos(1))) ' - ' char(g(pos(2)))]};
% txt={num2str(pos(2))}
end





