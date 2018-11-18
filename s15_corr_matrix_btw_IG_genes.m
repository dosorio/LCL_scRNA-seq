sortby="expr_level";
i_common_code;

a=gl123; a(7568)="TT";
i=find(extractBefore(a,3)=="IG");   % find all IG genes
g=gl123(i);

doimpute=false;
ethid="AFR";
switch ethid
    case "EUR"
        X=GM12878_expr; x=cellcycleGM12878; x1=i1; x2=i2; x3=i3;        
        GM12878_expr_G1=X(:,x1);  % G1[~,i]=sort(sum(A1,1),'descend'); A1=A1(:,i);
        GM12878_expr_G2M=X(:,x2); %[~,i]=sort(sum(A2,1),'descend'); A2=A2(:,i);
        GM12878_expr_S=X(:,x3);   %[~,i]=sort(sum(A3,1),'descend'); A3=A3(:,i);
        
        if doimpute
            M=i_Magic_impute(GM12878_expr_G1');
            M=M';
            A=M(i,:);
        else
            A=GM12878_expr_G1(i,:);
        end
        
    case "AFR"
        X=GM18502_expr; x=cellcycleGM18502; x1=j1; x2=j2; x3=j3;
        GM18502_expr_G1=X(:,x1);  % G1[~,i]=sort(sum(A1,1),'descend'); A1=A1(:,i);
        GM18502_expr_G2M=X(:,x2); % [~,i]=sort(sum(A2,1),'descend'); A2=A2(:,i);
        GM18502_expr_S=X(:,x3);   %[~,i]=sort(sum(A3,1),'descend'); A3=A3(:,i);
        
        if doimpute
            M=i_Magic_impute(GM18502_expr_G1');
            M=M';
            A=M(i,:);
        else
            A=GM18502_expr_G1(i,:);
        end
end

%%   

B=corr(A');
C=B>0.2;
C=C-diag(diag(C));

%%
i=sum(C)==0;
C2=C; C2(i,:)=[]; C2(:,i)=[];
g2=g; g2(i)=[];
bg1 = biograph(C2>0,g2);
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

function txt = myupdatefcn(~,event_obj,g)
% Customizes text of data tips
pos = event_obj.Position;
% idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
num2str(pos(2))
txt = {[char(g(pos(1))) ' - ' char(g(pos(2)))]};
% txt={num2str(pos(2))}
end





