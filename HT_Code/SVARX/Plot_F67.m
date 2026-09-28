%%% Nonlinear geopolitical shocks across three GPR series
%%% Based on the original first-step SVAR calculation


set(0,'defaultAxesFontName','Times');
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1.5);
rng('default');
rng(0);
addpath('Functions');
addpath('WorldGPR');

% Each script should write its data matrix to NEW.
dataScripts = {'NEWW2','NEWGPA2','NEWGPT2'};
shockNames = {'WORLD','ACT','THREAT'};
shockTags = {'W','A','T'};

for k = 1:numel(dataScripts)
    clear NEW
    run(dataScripts{k});

    vardata = NEW(:,1:3);
    pos_shock = 1;
    n = size(vardata,2);
    p = 2;
    c = 1;
    hor = 49;
    cum = 0;

    [D,~,~,~,~,~,~,~,~,~,eta] = chol_irf(vardata,n,p,c,hor);
    a = cumsum(D.^2,3);
    VVV = squeeze(a(:,pos_shock,:))./squeeze(sum(a,2));
    T1top = VVV(:,[1 13 25 49]);
    fprintf('\n%s shock\n',shockNames{k});
    MakeTable(T1top*100,1);

    finshock = eta(:,pos_shock);
    finshock_squared = eta(:,pos_shock).^2;
    T = length(finshock);
    dates = 1994 + (0:T-1)/4;

    figure('Name',sprintf('Geopolitical %s shocks',shockNames{k}), ...
           'NumberTitle','off');

    subplot(2,1,1);
    plot(dates,finshock,'-','LineWidth',2,'Color',[0.5 0.5 0.5]);
    axis tight;
    title(sprintf('Geopolitical %s shock, $u_{g%s,t}$', ...
          shockNames{k},shockTags{k}),'Interpreter','latex');
    set(gca,'FontSize',18);

    subplot(2,1,2);
    plot(dates,finshock_squared,'-','LineWidth',2,'Color',[0.2 0.2 0.2]);
    axis tight;
    title(sprintf('Geopolitical %s squared shock, $u_{g%s,t}^{2}$', ...
          shockNames{k},shockTags{k}),'Interpreter','latex');
    set(gca,'FontSize',18);
end
