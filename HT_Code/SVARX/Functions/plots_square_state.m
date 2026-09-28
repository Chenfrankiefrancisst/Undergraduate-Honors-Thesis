%% Plots by negative vs positive:

Shocknames={'Financial Shock - 1 standard deviation';'Financial Shock - 2 standard deviations'; 'Financial Shock - 4 standard deviations'};

magnitude=1;

C_state=magnitude.*C_wold(:,1:m:end)+(magnitude.^2).*C_wold(:,2:m:end)+magnitude.*C_wold(:,3:m:end);

[~,~,HighC_state,LowC_state]=bootstrapVARX_state(finaldata,news,hor,1,trend,iter,conf,p,q,n,magnitude,cum); % function that performs bootstrap
[~,~,HighC90_state,LowC90_state]=bootstrapVARX_state(finaldata,news,hor,1,trend,iter,90,p,q,n,magnitude,cum); % function that performs bootstrap

C_nostate=magnitude.*C_wold(:,1:m:end)-(magnitude.^2).*C_wold(:,2:m:end)+magnitude.*C_wold(:,3:m:end);

% [~,~,HighC_nostate,LowC_nostate]=bootstrapVARX_nostate(finaldata,news,hor,1,trend,iter,conf,p,q,n,magnitude,cum); % function that performs bootstrap
% [~,~,HighC90_nostate,LowC90_nostate]=bootstrapVARX_nostate(finaldata,news,hor,1,trend,iter,90,p,q,n,magnitude,cum); % function that performs bootstrap

set(0,'defaultAxesFontName', 'Times');
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1.5)

figure;
    
for k=[1:n;1:3:3*n] % only first shock
    
    subplot(n,3,k(2))
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

legend([s1,s2,p1,p2],{'90% confidence bands','68% confidence bands','IRF positive','IRF negative'},'FontSize',18,'Orientation','Horizontal')
    

magnitude=2;

C_state=magnitude.*C_wold(:,1:m:end)+(magnitude.^2).*C_wold(:,2:m:end)+magnitude.*C_wold(:,3:m:end);

[~,~,HighC_state,LowC_state]=bootstrapVARX_state(finaldata,news,hor,1,trend,iter,conf,p,q,n,magnitude,cum); % function that performs bootstrap
[~,~,HighC90_state,LowC90_state]=bootstrapVARX_state(finaldata,news,hor,1,trend,iter,90,p,q,n,magnitude,cum); % function that performs bootstrap

C_nostate=magnitude.*C_wold(:,1:m:end)-(magnitude.^2).*C_wold(:,2:m:end)+magnitude.*C_wold(:,3:m:end);

% [~,~,HighC_nostate,LowC_nostate]=bootstrapVARX_nostate(finaldata,news,hor,1,trend,iter,conf,p,q,n,magnitude,cum); % function that performs bootstrap
% [~,~,HighC90_nostate,LowC90_nostate]=bootstrapVARX_nostate(finaldata,news,hor,1,trend,iter,90,p,q,n,magnitude,cum); % function that performs bootstrap

for k=[1:n;1:3:3*n] % only first shock
    
    subplot(n,3,1+k(2))
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
        title(strcat(Shocknames{2}), 'FontSize', 18);
    end
    axis tight
    xlim([0 hor-1]);

end

magnitude=4;

C_state=magnitude.*C_wold(:,1:m:end)+(magnitude.^2).*C_wold(:,2:m:end)+magnitude.*C_wold(:,3:m:end);

[~,~,HighC_state,LowC_state]=bootstrapVARX_state(finaldata,news,hor,1,trend,iter,conf,p,q,n,magnitude,cum); % function that performs bootstrap
[~,~,HighC90_state,LowC90_state]=bootstrapVARX_state(finaldata,news,hor,1,trend,iter,90,p,q,n,magnitude,cum); % function that performs bootstrap

C_nostate=magnitude.*C_wold(:,1:m:end)-(magnitude.^2).*C_wold(:,2:m:end)+magnitude.*C_wold(:,3:m:end);

% [~,~,HighC_nostate,LowC_nostate]=bootstrapVARX_nostate(finaldata,news,hor,1,trend,iter,conf,p,q,n,magnitude,cum); % function that performs bootstrap
% [~,~,HighC90_nostate,LowC90_nostate]=bootstrapVARX_nostate(finaldata,news,hor,1,trend,iter,90,p,q,n,magnitude,cum); % function that performs bootstrap

for k=[1:n;1:3:3*n] % only first shock
    
    subplot(n,3,2+k(2))
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
        title(strcat(Shocknames{3}), 'FontSize', 18);
    end
    axis tight
    xlim([0 hor-1]);

end

