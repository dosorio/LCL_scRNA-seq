sortby="expr_level";
i_common_code;

ethid="EUR";
switch ethid
    case "EUR"
X=GM12878_expr; x=cellcycleGM12878; x1=i1; x2=i2; x3=i3;
    case "AFR"
X=GM18502_expr; x=cellcycleGM18502; x1=j1; x2=j2; x3=j3;
end


% a=gl123; a(7568)="TT";
% i=find(extractBefore(a,3)=="IG");   % find all IG genes
% X(i,:)=[];
% gl123(i)=[];
% gl123desc(i)=[];



%%
% figure;
% hold on
% for k=0.9:-0.05:0.5
%     plot(k, sum((sum(A>0,2)/size(A,2))>k),'ko');
% end

%%
%{
[xs,i]=sort(x);
figure;
imagesc(log10(1+X(1:2000,i)))
vline(find(xs=="G2M", 1 ),'y-')
vline(find(xs=="S", 1 ),'g-')
title('G1  |  G2/M  |   S')
colorbar;
xlabel('Cell')
ylabel('Gene')
%}

%%
figure;
A1=X(:,x1); [~,i]=sort(sum(A1,1),'descend'); A1=A1(:,i);
A2=X(:,x2); [~,i]=sort(sum(A2,1),'descend'); A2=A2(:,i);
A3=X(:,x3); [~,i]=sort(sum(A3,1),'descend'); A3=A3(:,i);

XA=[A1 A3 A2];
imagesc(log10(1+XA(1:200,:)))
% vline(find(x=="G2M", 1 ),'y-')
% vline(find(x=="S", 1 ),'y-')
vline(size(A1,2),'y-');
vline(size(A1,2)+size(A3,2),'y-');
title('G1  |  S  |   G2/M')
colorbar;
xlabel('Cell')
ylabel('Gene')
dt = datacursormode;
dt.UpdateFcn = {@myupdatefcn,gl123,gl123desc,0};


%%

% ethid="EUR"
switch ethid
    case "EUR"
        [~,i]=sort(sum(GMmix_eur_G1,1),'descend'); A1=GMmix_eur_G1(:,i);
        [~,i]=sort(sum(GMmix_eur_G2M,1),'descend'); A2=GMmix_eur_G2M(:,i);
        [~,i]=sort(sum(GMmix_eur_S,1),'descend'); A3=GMmix_eur_S(:,i);
    case "AFR"
        [~,i]=sort(sum(GMmix_afr_G1,1),'descend'); A1=GMmix_afr_G1(:,i);
        [~,i]=sort(sum(GMmix_afr_G2M,1),'descend'); A2=GMmix_afr_G2M(:,i);
        [~,i]=sort(sum(GMmix_afr_S,1),'descend'); A3=GMmix_afr_S(:,i);
end

figure;
XA=[A1 A3 A2];
imagesc(log10(1+XA(1:200,:)))
vline(size(A1,2),'y-');
vline(size(A1,2)+size(A3,2),'y-');
title('G1  |  S  |   G2/M')
colorbar;
xlabel('Cell')
ylabel('Gene')
dt = datacursormode;
dt.UpdateFcn = {@myupdatefcn,gl123,gl123desc,0};



%%
figure; 
subplot(2,1,1)
i=gl123=='RPS11';
stem(GMmix_eur_G1(i,:),'marker','none');
title([gl123(i) 'Ribosomal protein S11'])
ylim([0 700])
Entropy(GMmix_eur_G1(i,:)')
mean(GMmix_eur_G1(i,:))

subplot(2,1,2)
i=gl123=='IGHG1';
stem(GMmix_eur_G1(i,:),'marker','none');
title([gl123(i) 'Immunoglobulin Heavy Constant Gamma 1'])
ylim([0 700])
Entropy(GMmix_eur_G1(i,:)')
mean(GMmix_eur_G1(i,:))
%%
figure;
i=gl123=='MS4A1';
stem(GMmix_eur_G1(i,:),'marker','none');
title([gl123(i) gl123desc(i)])
Entropy(GMmix_eur_G1(i,:)')
web('https://dice-database.org/cells/B_CELL_NAIVE')

%%
function txt = myupdatefcn(~,event_obj,g,gd,id)
% Customizes text of data tips
pos = event_obj.Position;
% idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
num2str(pos(2)+id)
txt = {[num2str(pos(2)+id) ' - ' char(g(pos(2)+id)) ' - ' char(gd(pos(2)+id))]};
% txt={num2str(pos(2))}
end