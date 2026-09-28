figure;

for k=[1:n;1:n+m:n*n+n*m]
    
    subplot(4,2,k(1))
    b1=area(0:hor-1,[sum(vardec(:,k(2):k(2)+n+m-3),2) vardec(:,k(2)+n+m-2:k(2)+n+m-1)],'LineWidth',1); hold on; 
    set(gca,'FontSize',16)
    set(b1(1),'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
    set(b1(2),'FaceColor',[0.5 0.5 0.5],'EdgeColor','none')
    set(b1(3),'FaceColor',[0.2 0.2 0.2],'EdgeColor','none')
    title(strcat(VARnames{k(1)}), 'FontSize', 18,'Interpreter','Latex');
    xlim([0 hor-1]);
    ylim([0 1])
    
end

legend(b1,'Residual','$u_{ft}$','$g(u_{ft})$','Interpreter','Latex','Orientation','Horizontal','FontSize',22)
set(gca,'FontSize',16)
