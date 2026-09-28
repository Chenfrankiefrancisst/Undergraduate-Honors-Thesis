%% Plots by negative vs positive:

Shocknames={'Overall Financial Shock'};

magnitude=1;

C_pos=magnitude.*C_wold(:,1:m:end)+(magnitude).*C_wold(:,2:m:end);

[~,~,HighC_pos,LowC_pos]=bootstrapSVARX_pos(finaldata,news,hor,1,trend,iter,conf,p,q,n,magnitude,cum,pos_shock,restr,nonlin); % function that performs bootstrap
[~,~,HighC90_pos,LowC90_pos]=bootstrapSVARX_pos(finaldata,news,hor,1,trend,iter,90,p,q,n,magnitude,cum,pos_shock,restr,nonlin); % function that performs bootstrap

C_neg=magnitude.*C_wold(:,1:m:end)-(magnitude).*C_wold(:,2:m:end);

[~,~,HighC_neg,LowC_neg]=bootstrapSVARX_neg(finaldata,news,hor,1,trend,iter,conf,p,q,n,magnitude,cum,pos_shock,restr,nonlin); % function that performs bootstrap
[~,~,HighC90_neg,LowC90_neg]=bootstrapSVARX_neg(finaldata,news,hor,1,trend,iter,90,p,q,n,magnitude,cum,pos_shock,restr,nonlin); % function that performs bootstrap

set(0,'defaultAxesFontName', 'Times');
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1.5)
    
for k=[1:n;1:n] % only first shock
    
    subplot(n,1,k(2))
    s1=fill([0:hor-1 fliplr(0:hor-1)]' ,[HighC90_pos(:,k(1)); flipud(LowC90_pos(:,k(1)))],...
        colorBNDS,'EdgeColor','None'); hold on; 
    s2=fill([0:hor-1 fliplr(0:hor-1)]' ,[HighC_pos(:,k(1)); flipud(LowC_pos(:,k(1)))],...
        colorBNDS68,'EdgeColor','None'); hold on; %grid on; %grid minor;
    p1=plot(0:hor-1,C_pos(:,k(1)),'LineWidth',3.5,'Color','k'); hold on;
    p2=plot(0:hor-1,C_neg(:,k(1)),'LineWidth',3.5,'Color','r','LineStyle',':'); hold on;
    set(gca,'FontSize',18)
    line(get(gca,'Xlim'),[0 0],'Color','k','LineStyle','--','LineWidth',1.5); hold off;
    ylabel(strcat(VARnames{k(1)}), 'FontSize', 18);
    if k == 1
        title(strcat(Shocknames{1}), 'FontSize', 18);
    end
    axis tight
    xlim([0 hor-1]);

end

legend([s1,s2,p1,p2],{'90% confidence bands','68% confidence bands','IRF positive','IRF negative'},'FontSize',18,'Orientation','Horizontal')
    