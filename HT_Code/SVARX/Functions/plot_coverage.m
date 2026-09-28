figure;


for i=1:n*m
        subplot(3,2,i)
        p1=plot(0:hor-1,mean(coveragerate(:,i,:),3),'k','LineWidth',1.5,'LineStyle','--'); hold on;
        p2=plot(0:hor-1,mean(coveragerate_10000(:,i,:),3),'r','LineWidth',1.5,'LineStyle','-.'); hold on;
        z = plot(0:hor-1,zeros(hor,1)+0.68);
        set(z(1),'Color','b','LineWidth',1,'LineStyle','-');
        set(gca,'FontSize',12)
        axis([0 19 0 1])
        if i==1
        title(shock_names{i},'Interpreter','Latex')    
        end
        if i==2
        title(shock_names{i},'Interpreter','Latex')    
        end
end

legend([p1,p2,z],{'T=518','T=10000','68% coverage'},'FontSize',14,'Orientation','Horizontal')
 