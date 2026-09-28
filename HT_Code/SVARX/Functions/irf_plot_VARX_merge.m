function irf_plot_VARX_merge(MiddleD,HighD,LowD,HighD68,LowD68,VARnames,Shocknames,colorBNDS,colorBNDS68,hor,n,m,pos_shock)

set(0,'defaultAxesFontName', 'Times');
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1.5)

for i=1:m
    
for k=[1:m:n*m;1:n;pos_shock:n:n*n]
    
    subplot(n,m,k(1)+(i-1))
    s1=fill([0:hor-1 fliplr(0:hor-1)]' ,[HighD(:,k(1)+(i-1)); flipud(LowD(:,k(1)+(i-1)))],...
        colorBNDS,'EdgeColor','None'); hold on; 
    s2=fill([0:hor-1 fliplr(0:hor-1)]' ,[HighD68(:,k(1)+(i-1)); flipud(LowD68(:,k(1)+(i-1)))],...
        colorBNDS68,'EdgeColor','None'); hold on; %grid on; %grid minor;
    p1=plot(0:hor-1,MiddleD(:,k(1)+(i-1)),'LineWidth',3.5,'Color','k'); hold on;
    set(gca,'FontSize',18)
    line(get(gca,'Xlim'),[0 0],'Color','k','LineStyle','--','LineWidth',1.5); hold off;
    ylabel(strcat(VARnames{k(2)}), 'FontSize', 18);
    if k(1) == 1
        title(strcat(Shocknames{i}), 'FontSize', 18);
    end
    axis tight
    xlim([0 hor-1]);
    
end

end

legend([s1,s2,p1],{'90% confidence bands','68% confidence bands','IRF'},'FontSize',18,'Orientation','Horizontal')
    
end

