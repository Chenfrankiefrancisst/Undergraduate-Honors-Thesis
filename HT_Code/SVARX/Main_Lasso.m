
close all; clear; clc;

%% Settings and data upload:

set(0,'defaultAxesFontName', 'Times', ...
      'defaultTextInterpreter','latex', ...
      'defaultAxesTickLabelInterpreter','latex', ...
      'defaultLegendInterpreter','latex');
set(0,'defaultAxesLineStyleOrder','-|--|:', ...
      'defaultLineLineWidth',1.5)

rng('default');
rng(0);
addpath('Functions');
addpath('Data');

% Load dataset (must create NEW in workspace)
MDGPTW; 

vardata = NEW(:,1:9);  % T x 9
Yall    = vardata;
[T, n]  = size(vardata);

VARnames = { ...
    'GPR Threat'; ...
    'CPI';...
    'Policy Rate'; ...
    '1 Year Government Yield'; ...
    'Real Economic Activity'; ...
    'Oil Price'; ...
    'Gold Spot'; ...
    'Durables'; ...
    'SP 500' ...
};


pos_shock = 1;   % position of GPR shock in VAR
p_chol    = 1;   % # lags
p         = p_chol;  % for chol_irf / SVARX
c         = 1;   % intercept
hor       = 49;  % # horizons of IRFs (rows)
cum       = 0;   % IRFs not cumulated

% First-step SVAR to identify shocks
[D,~,~,~,~,~,~,~,~,~,eta] = chol_irf(vardata, n, p, c, hor);  % eta = Cholesky shocks

a     = cumsum(D.^2,3);
VVV   = squeeze(a(:,pos_shock,:))./squeeze(sum(a,2));
T1top = VVV(:,[1 13 25 49]);
MakeTable(T1top*100,1);

% Cholesky-identified GPR shocks (first structural shock)
finshock         = eta(:,pos_shock);
finshock_squared = finshock.^2;

% Plot Cholesky shocks
T_shock = length(finshock);
dates   = 1994 + (0:T_shock-1)/4;   % quarterly from 1994Q1

figure;
subplot(2,1,1)
plot(dates,finshock,'-','LineWidth',2,'Color',[0.5 0.5 0.5]); axis tight;
title('Cholesky GPR shock, $u_{GPR,t}$')
set(gca,'FontSize',16)
subplot(2,1,2)
plot(dates,finshock_squared,'-','LineWidth',2,'Color',[0.2 0.2 0.2]); axis tight;
title('Cholesky squared GPR shock, $u_{GPR,t}^2$')
set(gca,'FontSize',16)


finaldata = vardata(p+1:end,:);     % endogenous variables of VARX
news      = [finshock, finshock_squared];  % exogenous (linear, squared)

nonlin = 1;   % square
restr  = 2;   % same as your code
m      = size(news,2);  % # exogenous vars
q      = 0;   % # lags of exogenous

SVARX;  

% Standard Forni-style IRF plot (Cholesky only)
colorBNDS   = [0.9 0.9 0.9];
colorBNDS68 = [0.68 0.68 0.68];

Shocknames = {'GPR Shock'; 'GPR Shock Squared'};

irf_plot_VARX(C_wold,HighC90,LowC90,HighC,LowC, ...
              VARnames,Shocknames,colorBNDS,colorBNDS68,hor,n,m,pos_shock);


s_chol   = 2 * std(finshock);   % size of identified shock
H        = hor - 1;        % horizon index from 0..H
horizons = (0:H)';

m_exo = 2;   % linear, squared

% Storage for Cholesky IRFs and bands
chol_lin_irf   = zeros(hor,n);
chol_lin_hi90  = zeros(hor,n);
chol_lin_lo90  = zeros(hor,n);
chol_lin_hi68  = zeros(hor,n);
chol_lin_lo68  = zeros(hor,n);

chol_sq_irf    = zeros(hor,n);
chol_sq_hi90   = zeros(hor,n);
chol_sq_lo90   = zeros(hor,n);
chol_sq_hi68   = zeros(hor,n);
chol_sq_lo68   = zeros(hor,n);

chol_tot_irf   = zeros(hor,n);
chol_tot_hi90  = zeros(hor,n);
chol_tot_lo90  = zeros(hor,n);
chol_tot_hi68  = zeros(hor,n);
chol_tot_lo68  = zeros(hor,n);

for v_idx = 1:n
    col_lin = (v_idx-1)*m_exo + 1;
    col_sq  = (v_idx-1)*m_exo + 2;

    % ---- linear channel (scaled) ----
    lin_irf   = s_chol  * C_wold(:, col_lin);
    lin_hi90  = s_chol  * HighC90(:, col_lin);
    lin_lo90  = s_chol  * LowC90(:,  col_lin);
    lin_hi68  = s_chol  * HighC(:,    col_lin);
    lin_lo68  = s_chol  * LowC(:,     col_lin);

    % ---- squared channel (scaled) ----
    sq_irf    = (s_chol^2) * C_wold(:, col_sq);
    sq_hi90   = (s_chol^2) * HighC90(:, col_sq);
    sq_lo90   = (s_chol^2) * LowC90(:,  col_sq);
    sq_hi68   = (s_chol^2) * HighC(:,    col_sq);
    sq_lo68   = (s_chol^2) * LowC(:,     col_sq);

    % Store linear
    chol_lin_irf(:,v_idx)  = lin_irf;
    chol_lin_hi90(:,v_idx) = lin_hi90;
    chol_lin_lo90(:,v_idx) = lin_lo90;
    chol_lin_hi68(:,v_idx) = lin_hi68;
    chol_lin_lo68(:,v_idx) = lin_lo68;

    % Store squared
    chol_sq_irf(:,v_idx)   = sq_irf;
    chol_sq_hi90(:,v_idx)  = sq_hi90;
    chol_sq_lo90(:,v_idx)  = sq_lo90;
    chol_sq_hi68(:,v_idx)  = sq_hi68;
    chol_sq_lo68(:,v_idx)  = sq_lo68;

    % ---- total effect: combine lin & sq bands in quadrature ----
    tot_irf = lin_irf + sq_irf;

    % 90% band widths
    lin_up90 = lin_hi90 - lin_irf;
    lin_dn90 = lin_irf  - lin_lo90;
    sq_up90  = sq_hi90  - sq_irf;
    sq_dn90  = sq_irf   - sq_lo90;

    tot_up90 = sqrt(lin_up90.^2 + sq_up90.^2);
    tot_dn90 = sqrt(lin_dn90.^2 + sq_dn90.^2);

    tot_hi90 = tot_irf + tot_up90;
    tot_lo90 = tot_irf - tot_dn90;

    % 68% band widths
    lin_up68 = lin_hi68 - lin_irf;
    lin_dn68 = lin_irf  - lin_lo68;
    sq_up68  = sq_hi68  - sq_irf;
    sq_dn68  = sq_irf   - sq_lo68;

    tot_up68 = sqrt(lin_up68.^2 + sq_up68.^2);
    tot_dn68 = sqrt(lin_dn68.^2 + sq_dn68.^2);

    tot_hi68 = tot_irf + tot_up68;
    tot_lo68 = tot_irf - tot_dn68;

    chol_tot_irf(:,v_idx)   = tot_irf;
    chol_tot_hi90(:,v_idx)  = tot_hi90;
    chol_tot_lo90(:,v_idx)  = tot_lo90;
    chol_tot_hi68(:,v_idx)  = tot_hi68;
    chol_tot_lo68(:,v_idx)  = tot_lo68;
end

p_lasso = 1;


t_start = p_lasso + 1;   % earliest t where lags exist

gpr_col = 1;
Y_dep = Yall(t_start:T, gpr_col);   

X_reg = [];
for lag = 1:p_lasso
    X_reg = [X_reg, Yall(t_start-lag:T-lag, :)];  
end

[X_std, muX, sigX] = zscore(X_reg);
[Y_std, muY, sigY] = zscore(Y_dep);

[B_lasso, FitInfo] = lasso(X_std, Y_std, 'CV', 10);
idx        = FitInfo.IndexMinMSE;
beta_std   = B_lasso(:, idx);
alpha_std  = FitInfo.Intercept(idx);
Yhat_std   = X_std * beta_std + alpha_std;
Yhat       = Yhat_std * sigY + muY;

u_lasso = Y_dep - Yhat;      

fprintf('LASSO shock: var = %.6f, std = %.6f\n', var(u_lasso), std(u_lasso));

shock_lin = u_lasso;
shock_sq  = u_lasso.^2;

Y_varx = Yall(t_start:T, :);      
T_eff  = size(Y_varx,1);

const = ones(T_eff,1);

X_lags = [];
for lag = 1:p_lasso
    X_lags = [X_lags, Yall(t_start-lag:T-lag, :)];
end

X_shock = [shock_lin, shock_sq];

X_varx      = [const, X_lags, X_shock];    
Y_dep_varx  = Y_varx;                     

B = X_varx \ Y_dep_varx;                 

k      = size(X_varx,2);
k_lags = 1 + n*p_lasso;

A = cell(p_lasso,1);
row_start = 2;
for lag = 1:p_lasso
    rows = row_start : row_start + n - 1;
    A{lag} = B(rows,:).';         % n x n
    row_start = row_start + n;
end

idx_lin_lasso = k_lags + 1;
idx_sq_lasso  = k_lags + 2;
B_lin_lasso   = B(idx_lin_lasso,:).';   % n x 1
B_sq_lasso    = B(idx_sq_lasso,:).';    % n x 1

fprintf('LASSO VARX shock coeff (eq 1): lin = %.4f, sq = %.4f\n', ...
        B_lin_lasso(1), B_sq_lasso(1));


np_lasso = n * p_lasso;
F_lasso  = zeros(np_lasso, np_lasso);

for lag = 1:p_lasso
    F_lasso(1:n, (lag-1)*n+1:lag*n) = A{lag};
end
if p_lasso > 1
    F_lasso(n+1:end,1:(p_lasso-1)*n) = eye((p_lasso-1)*n);
end

G_lasso = zeros(np_lasso, 2);
G_lasso(1:n,1) = B_lin_lasso;
G_lasso(1:n,2) = B_sq_lasso;

H        = hor - 1;    
horizons = (0:H)';

s_lasso = 2 * std(shock_lin);

IRF_lin_L   = zeros(H+1,n);
IRF_sq_L    = zeros(H+1,n);
IRF_tot_L   = zeros(H+1,n);

S_lin0 = G_lasso(:,1)*s_lasso;
S_sq0  = G_lasso(:,2)*(s_lasso^2);
S_tot0 = S_lin0 + S_sq0;

IRF_lin_L(1,:) = S_lin0(1:n).';
IRF_sq_L(1,:)  = S_sq0(1:n).';
IRF_tot_L(1,:) = S_tot0(1:n).';

S_lin = S_lin0;
S_sq  = S_sq0;
S_tot = S_tot0;

for h = 1:H
    S_lin = F_lasso * S_lin;
    S_sq  = F_lasso * S_sq;
    S_tot = F_lasso * S_tot;

    IRF_lin_L(h+1,:) = S_lin(1:n).';
    IRF_sq_L(h+1,:)  = S_sq(1:n).';
    IRF_tot_L(h+1,:) = S_tot(1:n).';
end


Bboots = 500;
IRF_lin_boot = zeros(Bboots, H+1, n);
IRF_sq_boot  = zeros(Bboots, H+1, n);
IRF_tot_boot = zeros(Bboots, H+1, n);

% VARX residuals
resid = Y_dep_varx - X_varx * B;
resid = resid - mean(resid);

fprintf('Bootstrapping LASSO VARX for IRF bands (%d draws)...\n', Bboots);

for b = 1:Bboots
    idx = randi(T_eff,T_eff,1);
    e_star = resid(idx,:);

    Y_star = X_varx * B + e_star;

    B_star = X_varx \ Y_star;

    row_start = 2;
    A_star = cell(p_lasso,1);
    for lag = 1:p_lasso
        rows = row_start : row_start + n - 1;
        A_star{lag} = B_star(rows,:).';
        row_start = row_start + n;
    end

    B_lin_star = B_star(k_lags+1,:).';
    B_sq_star  = B_star(k_lags+2,:).';

    F_star = zeros(np_lasso,np_lasso);
    for lag = 1:p_lasso
        F_star(1:n, (lag-1)*n+1:lag*n) = A_star{lag};
    end
    if p_lasso > 1
        F_star(n+1:end,1:(p_lasso-1)*n) = eye((p_lasso-1)*n);
    end

    G_star = zeros(np_lasso,2);
    G_star(1:n,1) = B_lin_star;
    G_star(1:n,2) = B_sq_star;

    S_lin_b = G_star(:,1)*s_lasso;
    S_sq_b  = G_star(:,2)*(s_lasso^2);
    S_tot_b = S_lin_b + S_sq_b;

    lin_tmp = zeros(H+1,n);
    sq_tmp  = zeros(H+1,n);
    tot_tmp = zeros(H+1,n);

    lin_tmp(1,:) = S_lin_b(1:n).';
    sq_tmp(1,:)  = S_sq_b(1:n).';
    tot_tmp(1,:) = S_tot_b(1:n).';

    for h = 1:H
        S_lin_b = F_star * S_lin_b;
        S_sq_b  = F_star * S_sq_b;
        S_tot_b = F_star * S_tot_b;

        lin_tmp(h+1,:) = S_lin_b(1:n).';
        sq_tmp(h+1,:)  = S_sq_b(1:n).';
        tot_tmp(h+1,:) = S_tot_b(1:n).';
    end

    IRF_lin_boot(b,:,:) = lin_tmp;
    IRF_sq_boot(b,:,:)  = sq_tmp;
    IRF_tot_boot(b,:,:) = tot_tmp;
end

fprintf('LASSO bootstrap completed.\n');

% LASSO medians and CIs
lasso_lin_med  = squeeze(median(IRF_lin_boot,1));
lasso_lin_hi68 = squeeze(prctile(IRF_lin_boot,84,1));
lasso_lin_lo68 = squeeze(prctile(IRF_lin_boot,16,1));
lasso_lin_hi90 = squeeze(prctile(IRF_lin_boot,95,1));
lasso_lin_lo90 = squeeze(prctile(IRF_lin_boot,5, 1));

lasso_sq_med   = squeeze(median(IRF_sq_boot,1));
lasso_sq_hi68  = squeeze(prctile(IRF_sq_boot,84,1));
lasso_sq_lo68  = squeeze(prctile(IRF_sq_boot,16,1));
lasso_sq_hi90  = squeeze(prctile(IRF_sq_boot,95,1));
lasso_sq_lo90  = squeeze(prctile(IRF_sq_boot,5, 1));

lasso_tot_med  = squeeze(median(IRF_tot_boot,1));
lasso_tot_hi68 = squeeze(prctile(IRF_tot_boot,84,1));
lasso_tot_lo68 = squeeze(prctile(IRF_tot_boot,16,1));
lasso_tot_hi90 = squeeze(prctile(IRF_tot_boot,95,1));
lasso_tot_lo90 = squeeze(prctile(IRF_tot_boot,5, 1));


% Colors
chol90_col = [0.85 0.85 0.85];
chol68_col = [0.65 0.65 0.65];
lass90_col = [0.80 0.90 1.00];
lass68_col = [0.60 0.80 1.00];

figure;
for i = 1:n
    subplot(3,3,i); hold on; box on;
    set(gca,'FontSize',12);

    fill([horizons; flipud(horizons)], ...
         [chol_lin_hi90(:,i); flipud(chol_lin_lo90(:,i))], ...
         chol90_col, 'EdgeColor','none');
    fill([horizons; flipud(horizons)], ...
         [chol_lin_hi68(:,i); flipud(chol_lin_lo68(:,i))], ...
         chol68_col, 'EdgeColor','none');

    fill([horizons; flipud(horizons)], ...
         [lasso_lin_hi90(:,i); flipud(lasso_lin_lo90(:,i))], ...
         lass90_col, 'EdgeColor','none', 'FaceAlpha',0.7);
    fill([horizons; flipud(horizons)], ...
         [lasso_lin_hi68(:,i); flipud(lasso_lin_lo68(:,i))], ...
         lass68_col, 'EdgeColor','none', 'FaceAlpha',0.7);

    plot(horizons, chol_lin_irf(:,i), 'k-', 'LineWidth',1.8);
    plot(horizons, lasso_lin_med(:,i), 'b--', 'LineWidth',1.8);

    yline(0,'k-');

    title(VARnames{i});

    if i == 1
        legend({'Chol 90\%','Chol 68\%','LASSO 90\%','LASSO 68\%', ...
                'Chol IRF','LASSO IRF'}, ...
                'Location','best');
    end
end
sgtitle('Linear GPR Shock: Cholesky (black/gray) vs LASSO (blue)');

figure;
for i = 1:n
    subplot(3,3,i); hold on; box on;
    set(gca,'FontSize',12);

    fill([horizons; flipud(horizons)], ...
         [chol_sq_hi90(:,i); flipud(chol_sq_lo90(:,i))], ...
         chol90_col,'EdgeColor','none');
    fill([horizons; flipud(horizons)], ...
         [chol_sq_hi68(:,i); flipud(chol_sq_lo68(:,i))], ...
         chol68_col,'EdgeColor','none');

    fill([horizons; flipud(horizons)], ...
         [lasso_sq_hi90(:,i); flipud(lasso_sq_lo90(:,i))], ...
         lass90_col,'EdgeColor','none','FaceAlpha',0.7);
    fill([horizons; flipud(horizons)], ...
         [lasso_sq_hi68(:,i); flipud(lasso_sq_lo68(:,i))], ...
         lass68_col,'EdgeColor','none','FaceAlpha',0.7);

    plot(horizons, chol_sq_irf(:,i), 'k-', 'LineWidth',1.8);
    plot(horizons, lasso_sq_med(:,i), 'b--', 'LineWidth',1.8);

    yline(0,'k-');
    title(VARnames{i});
end
sgtitle('Squared GPR Shock: Cholesky (black/gray) vs LASSO (blue)');

figure;
for i = 1:n
    subplot(3,3,i); hold on; box on;
    set(gca,'FontSize',12);

    fill([horizons; flipud(horizons)], ...
         [chol_tot_hi90(:,i); flipud(chol_tot_lo90(:,i))], ...
         chol90_col, 'EdgeColor','none');
    fill([horizons; flipud(horizons)], ...
         [chol_tot_hi68(:,i); flipud(chol_tot_lo68(:,i))], ...
         chol68_col, 'EdgeColor','none');

    fill([horizons; flipud(horizons)], ...
         [lasso_tot_hi90(:,i); flipud(lasso_tot_lo90(:,i))], ...
         lass90_col, 'EdgeColor','none', 'FaceAlpha',0.7);
    fill([horizons; flipud(horizons)], ...
         [lasso_tot_hi68(:,i); flipud(lasso_tot_lo68(:,i))], ...
         lass68_col, 'EdgeColor','none', 'FaceAlpha',0.7);

    plot(horizons, chol_tot_irf(:,i), 'k-', 'LineWidth',1.8);
    plot(horizons, lasso_tot_med(:,i), 'b--', 'LineWidth',1.8);

    yline(0,'k-');
    title(VARnames{i});
end
sgtitle('Total Effect (Linear + Squared): Cholesky (black/gray) vs LASSO (blue)');
