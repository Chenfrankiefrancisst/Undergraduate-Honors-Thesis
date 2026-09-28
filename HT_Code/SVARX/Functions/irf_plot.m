function irf_plot(MiddleD,HighD,LowD,HighD68,LowD68,VARnames,Shocknames,colorBNDS,colorBNDS68,hor,n)

set(0,'defaultAxesFontName', 'Times');
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1.5)

figure;

for k=1:n % variable k
    for j=1:n % shock j
    subplot(n,n,j+n*k-n)
    fill([0:hor-1 fliplr(0:hor-1)]' ,[HighD(:,j+n*k-n); flipud(LowD(:,j+n*k-n))],...
        colorBNDS,'EdgeColor','None'); hold on;
    fill([0:hor-1 fliplr(0:hor-1)]' ,[HighD68(:,j+n*k-n); flipud(LowD68(:,j+n*k-n))],...
        colorBNDS68,'EdgeColor','None'); hold on; %grid on; %grid minor;
    %plot(0:hor-1,HighD(:,j+n*k-n),'LineWidth',1.5,'Color','k'); hold on;
    %plot(0:hor-1,HighD68(:,j+n*k-n),'LineWidth',1.5,'Color','k'); hold on;
    plot(0:hor-1,MiddleD(:,j+n*k-n),'LineWidth',3.5,'Color','k'); hold on;
    %plot(0:hor-1,LowD(:,j+n*k-n),'LineWidth',1.5,'Color','k'); hold on;
    %plot(0:hor-1,LowD68(:,j+n*k-n),'LineWidth',1.5,'Color','k'); hold on;
    set(gca,'FontSize',18)
    line(get(gca,'Xlim'),[0 0],'Color',[1 0 0],'LineStyle','-','LineWidth',1.5); hold off;
    
    if j ==1
        ylabel(strcat(VARnames{k}), 'FontSize', 18);
    end
    
    if k == 1
        title(strcat(Shocknames{j}), 'FontSize', 18);
    end
    xlim([0 hor-1]);
    end
end

end