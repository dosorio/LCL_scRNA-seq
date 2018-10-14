load res
close all

i=Td.ResultDiffDisp=="NoDiff" | Tm.ResultDiffMean=="NoDiff"; % | Td.ResultDiffDisp=="ExcludedFromTesting";
figure;
% plot(Tm.MeanLog2FC(i),Td.DispLog2FC(i),'o','color',[.5 .5 .5])
% hold on
% plot(Tm.MeanLog2FC(~i),Td.DispLog2FC(~i),'o')
plot(Tm.MeanLog2FC,Td.DispLog2FC,'o')
dt = datacursormode;
dt.UpdateFcn = {@myupdatefcn,Td.GeneName};



function txt = myupdatefcn(~,event_obj,g)
% Customizes text of data tips
%pos = event_obj.Position;
idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
% num2str(pos(2)+id)
% txt = {[char(g(pos(2)+id)) ' - ' char(gd(pos(2)+id))]};
txt={g(idx)};
end