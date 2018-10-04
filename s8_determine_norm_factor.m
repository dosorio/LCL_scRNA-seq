global GM12878_expr GM18502_expr GMmix_eur GMmix_afr
sortby="expr_level";
i_common_code;
close all
%%
ifunc=@mean;

median_eur=ifunc(GM12878_expr,2)*1000;
median_afr=ifunc(GM18502_expr,2)*1000;
median_mix=ifunc(GMmix_expr,2)*1000;
median_mix_eur=ifunc(GMmix_eur,2)*1000;
median_mix_afr=ifunc(GMmix_afr,2)*1000;

%%

figure;
subplot(3,3,1)
i_loglogs(median_afr,median_eur,'Afr','Eur',gl123);
subplot(3,3,2)
i_loglogs(median_mix_afr,median_mix_eur,'mixAfr','mixEur',gl123);
subplot(3,3,3)
i_loglogs(median_afr,median_mix_afr,'Afr','mixAfr',gl123);
subplot(3,3,4)
i_loglogs(median_eur,median_mix_eur,'Eur','mixEur',gl123);


function i_loglogs(m1,m2,lb1,lb2,g)
    loglog(m1,m2,'.');
    [r]=corr(m1,m2,'type','s');
    title(r); xlabel(lb1); ylabel(lb2)
    hold on
    loglog([0.1 1e5], [0.1 1e5], 'r')
    %loglog([1e3 1e6], [1e3 1e6], 'r')
    xlim([0.1 1e6])
    ylim([0.1 1e6])    
    dt = datacursormode;
    dt.UpdateFcn = {@myupdatefcn,g};
end

function i_loglogs2(m1,m2,lb1,lb2,g)
    loglog(m1,m2,'.');
    [r]=corr(m1,m2,'type','s');
    title(r); xlabel(lb1); ylabel(lb2)
    hold on
    %loglog([0.01 1e5], [0.01 1e5], 'r')
    loglog([1e3 1e6], [1e3 1e6], 'r')
    %xlim([0.01 1e6])
    %ylim([0.01 1e6])    
    dt = datacursormode;
    dt.UpdateFcn = {@myupdatefcn,g};
end


function txt = myupdatefcn(~,event_obj,g)
% Customizes text of data tips
% pos = event_obj.Position;
idx = event_obj.DataIndex;
i_plotsiglegene(idx,g);
txt = {g(idx)};

% txt = {['X: ',num2str(pos(1))],...
%        ['Y: ',num2str(pos(2))],...
%        ['I: ',num2str(idx)],...
%        ['G: ',char(g(idx))]};
end

function i_plotsiglegene(idx,g)
global GM12878_expr GM18502_expr GMmix_eur GMmix_afr
M1=GM12878_expr;
M2=GM18502_expr;
M3=GMmix_eur;
M4=GMmix_afr;

figure;
subplot(4,1,1)
stem((M1(idx,:)),'marker','none','color','k');
% xlim([0 length(M1(idx,:))])
xlim([0 5000])
title('GM12878 EUR')
y1=ylim;

subplot(4,1,2)
stem((M2(idx,:)),'marker','none','color','k');
xlim([0 5000])
%xlim([0 length(M2(idx,:))])
title('GM18502 AFR')
y2=ylim;

subplot(4,1,3)
stem((M3(idx,:)),'marker','none','color','k');
xlim([0 5000])
vline(size(M3,2),'r-')
%xlim([0 length(M3(idx,:))])
title('GMmix EUR')
ylim(y1);

subplot(4,1,4)
stem((M4(idx,:)),'marker','none','color','k');
vline(size(M4,2),'r-')
xlim([0 5000])
%xlim([0 length(M4(idx,:))])
title('GMmix AFR')
ylim(y2)
xlabel(sprintf('%s',g(idx)));
end