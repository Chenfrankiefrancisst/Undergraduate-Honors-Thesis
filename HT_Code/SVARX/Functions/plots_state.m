%% Plots by negative vs positive:

Shocknames={'Financial Shock'};

magnitude=1;

C_state=magnitude.*C_wold(:,1:m:end)+(magnitude.^2).*C_wold(:,2:m:end)+magnitude.*C_wold(:,3:m:end);

[~,~,HighC_state,LowC_state]=bootstrapSVARX_state(finaldata,news,hor,1,trend,iter,conf,p,q,n,magnitude,cum,nonlin); % function that performs bootstrap
[~,~,HighC90_state,LowC90_state]=bootstrapSVARX_state(finaldata,news,hor,1,trend,iter,90,p,q,n,magnitude,cum,nonlin); % function that performs bootstrap

C_nostate=magnitude.*C_wold(:,1:m:end)+(magnitude.^2).*C_wold(:,2:m:end);

set(0,'defaultAxesFontName', 'Times');
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1.5)

figure;
    
for k=[1:n;1:n] % only first shock
    
    subplot(n,1,k(2))
    s1=fill([0:hor-1 fliplr(0:hor-1)]' ,[HighC90_state(:,k(1)); flipud(LowC90_state(:,k(1)))],...
        colorBNDS,'EdgeColor','None'); hold on; 
    s2=fill([0:hor-1 fliplr(0:hor-1)]' ,[HighC_state(:,k(1)); flipud(LowC_state(:,k(1)))],...
        colorBNDS68,'EdgeColor','None'); hold on; %grid on; %grid minor;
    p1=plot(0:hor-1,C_state(:,k(1)),'LineWidth',3.5,'Color','k'); hold on;
    p2=plot(0:hor-1,C_nostate(:,k(1)),'LineWidth',3.5,'Color','r','LineStyle',':'); hold on;
    set(gca,'FontSize',18)
    line(get(gca,'Xlim'),[0 0],'Color','k','LineStyle','--','LineWidth',1.5); hold off;
    ylabel(strcat(VARnames{k(1)}), 'FontSize', 18);
    if k == 1
        title(strcat(Shocknames{1}), 'FontSize', 18);
    end
    axis tight
    xlim([0 hor-1]);

end

legend([s1,s2,p1,p2],{'90% confidence bands','68% confidence bands','IRF state','IRF no state'},'FontSize',18,'Orientation','Horizontal')
    

