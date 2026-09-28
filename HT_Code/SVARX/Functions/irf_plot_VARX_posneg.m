function irf_plot_VARX_posneg(MiddleD_1,HighD_1,LowD_1,HighD68_1,LowD68_1,MiddleD_2,HighD_2,LowD_2,HighD68_2,LowD68_2,VARnames,Shocknames,colorBNDS,colorBNDS68,hor,n,m)

set(0,'defaultAxesFontName', 'Times');
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1.5)

figure;
    
for k=[1:n;1:n] % only first shock
    
    subplot(n,1,k(1))
    s1=fill([0:hor-1 fliplr(0:hor-1)]' ,[HighD_1(:,k(1)); flipud(LowD_1(:,k(1)))],...
        colorBNDS,'EdgeColor','None'); hold on; 
    s2=fill([0:hor-1 fliplr(0:hor-1)]' ,[HighD68_1(:,k(1)); flipud(LowD68_1(:,k(1)))],...
        colorBNDS68,'EdgeColor','None'); hold on; %grid on; %grid minor;
    p1=plot(0:hor-1,MiddleD_1(:,k(1)),'LineWidth',3.5,'Color','k'); hold on;
%     fill([0:hor-1 fliplr(0:hor-1)]' ,[HighD_2(:,k(1)); flipud(LowD_2(:,k(1)))],...
%          colorBNDS,'EdgeColor','None'); hold on; 
%     fill([0:hor-1 fliplr(0:hor-1)]' ,[HighD68_2(:,k(1)); flipud(LowD68_2(:,k(1)))],...
%          colorBNDS68,'EdgeColor','None'); hold on; %grid on; %grid minor;
    p2=plot(0:hor-1,MiddleD_2(:,k(1)),'LineWidth',3.5,'Color','r','LineStyle',':'); hold on;
    set(gca,'FontSize',18)
    line(get(gca,'Xlim'),[0 0],'Color','k','LineStyle','--','LineWidth',1.5); hold off;
    ylabel(strcat(VARnames{k(2)}), 'FontSize', 18);
    if k == 1
        title(strcat(Shocknames{1}), 'FontSize', 18);
    end
    axis tight
    xlim([0 hor-1]);

end

legend([s1,s2,p1,p2],{'90% confidence bands','68% confidence bands','IRF negative','IRF positive'},'FontSize',18,'Orientation','Horizontal')
    
end
