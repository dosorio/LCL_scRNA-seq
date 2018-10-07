sortby="expr_level";
i_common_code;
close all
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
startgid=1;
viewgnum=200;

h=figure;
% XA=[A1 A3 A2 B1 B3 B2];
XA=[A1 B1];
imagesc(log10(1+XA(startgid:startgid+viewgnum,:)))
% vline(find(x=="G2M", 1 ),'y-')
% vline(find(x=="S", 1 ),'y-')
vline(size(A1,2),'y-');
% vline(size(A1,2)+size(A3,2),'y-');
% title('G1  |  S  |   G2/M')
yticks(1:viewgnum)
yticklabels(gl123(startgid:startgid+viewgnum))
%ytickangle(45)

title('EUR G1  |  AFR G1')
colorbar;
xlabel('Cell')
ylabel('Gene')
dt = datacursormode;
dt.UpdateFcn = {@myupdatefcn,gl123,gl123desc,startgid-1};


function txt = myupdatefcn(~,event_obj,g,gd,id)
% Customizes text of data tips
pos = event_obj.Position;
% idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
% num2str(pos(2)+id)
gidx=pos(2)+id;
txt = {[num2str(gidx) ' - ' char(g(gidx)) ' - ' char(gd(gidx))]};
% txt={num2str(pos(2))}
end