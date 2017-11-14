function ploterr(n,T,Ydata,MC_data,Yname,Ylabel,Linestyle)
figure('Position',[0,0,1200,400]);
axes1 = axes('YScale','log','YMinorTick','on');
set(gca,'FontSize',20);
set(gca,'LooseInset',get(gca,'TightInset'));
hold(axes1,'all');
for i = 1:n
    Ydata_err(i,:) = abs(Ydata(i,:)-MC_data);
end
semilogy1 = semilogy(T,Ydata_err,'Parent',axes1,'Color',[0 0 0]);
for i = 1:n
    set(semilogy1(i),'LineStyle',Linestyle{i},'DisplayName',Yname{i});
end
legend(axes1,'show','Location','NorthEast');
ylabel(Ylabel);
xlabel('t(s)');