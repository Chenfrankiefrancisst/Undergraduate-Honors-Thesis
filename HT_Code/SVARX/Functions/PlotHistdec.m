% Plot Historical:

dates=dates(p+1:end);

figure; 

for i=1:n
    
subplot(4,2,i)
b1=bar(dates,histdecpos(:,[i i+n i+2*n]),'stack');
set(b1(1),'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
set(b1(2),'FaceColor',[0.5 0.5 0.5],'EdgeColor','none')
set(b1(3),'FaceColor',[0.2 0.2 0.2],'EdgeColor','none')
hold on;
b2=bar(dates,histdecneg(:,[i i+n i+2*n]),'stack');
set(b2(1),'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
set(b2(2),'FaceColor',[0.5 0.5 0.5],'EdgeColor','none')
set(b2(3),'FaceColor',[0.2 0.2 0.2],'EdgeColor','none')
hold on;
h3=plot(dates,shock_contr(:,i),'-','LineWidth',2,'Color',[0.100 0.1 0.1]);hold on;
title(VARnames{i},'Interpreter','Latex')
set(gca,'FontSize',18)

end

legend(b1,'Residual','$u_{ft}$','$g(u_{ft})$','Interpreter','Latex','Orientation','Horizontal','FontSize',22)