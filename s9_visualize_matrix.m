sortby="expr_level";
i_common_code;

ethid="EUR";
switch ethid
    case "EUR"
X=GM12878_expr; x=cellcycleGM12878; x1=i1; x2=i2; x3=i3;
    case "AFR"
X=GM18502_expr; x=cellcycleGM18502; x1=j1; x2=j2; x3=j3;
end

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
imagesc(log10(1+XA(201:300,:)))
% vline(find(x=="G2M", 1 ),'y-')
% vline(find(x=="S", 1 ),'y-')
vline(size(A1,2),'y-');
vline(size(A1,2)+size(A3,2),'y-');
title('G1  |  S  |   G2/M')
colorbar;
xlabel('Cell')
ylabel('Gene')
dt = datacursormode;
dt.UpdateFcn = {@myupdatefcn,gl123,gl123desc,200};


function txt = myupdatefcn(~,event_obj,g,gd,id)
% Customizes text of data tips
pos = event_obj.Position;
% idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
num2str(pos(2)+id)
txt = {[char(g(pos(2)+id)) ' - ' char(gd(pos(2)+id))]};
% txt={num2str(pos(2))}
end