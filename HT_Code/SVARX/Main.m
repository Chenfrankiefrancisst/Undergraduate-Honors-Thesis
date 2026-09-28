% GPR shock comparison: linear, quadratic, and total IRFs

clear; clc; close all;
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1.5)
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
rng('default');
rng(0);
addpath('Functions');
addpath('Data');

% Aligned quarterly data
load('AlignedData_1994Q2_2023Q2.mat','NEW_28vars','NNEW');

% Same eight macro series for each shock
Macro = NNEW(:,2:9);
n = 9;
p = 1;
c = 1;
hor = 49;
cum = 0;
pos_shock = 1;
nonlin = 1;
restr = 2;
q = 0;
tvec = 0:hor-1;

% Column order in NEW_28vars
CountryNames = { ...
    'GPR(World)','GPR-Threat','GPR-Act', ...
    'Australia','Brazil','Canada','Switzerland','China','Germany','Egypt', ...
    'France','United Kingdom','India','Israel','Japan','South Korea', ...
    'Mexico','Netherlands','Russia','Saudi Arabia','Sweden','Turkey', ...
    'Taiwan','Ukraine','United States','Venezuela','Vietnam','South Africa' ...
};

% Compare these three GPR shocks
target_names = {'GPR(World)','GPR-Act','GPR-Threat'};
target_colors = [0 0 0;
                 1 0 0;
                 0 0 1];
numT = numel(target_names);

% Panels in each figure
VARnames = { ...
    'GPR World'; ...
    'CPIU'; ...
    'Shadow Rate'; ...
    '1Y Gov Yield'; ...
    'Gold Demand';...
    'Gold Production World';...
    'Real Activity';...
    'Gold Spot';...
    'Durables';...
    };
LinIRF = zeros(hor,n,numT);
SqIRF = zeros(hor,n,numT);
TotIRF = zeros(hor,n,numT);
LinHi68 = zeros(hor,n,numT);
LinLo68 = zeros(hor,n,numT);
LinHi90 = zeros(hor,n,numT);
LinLo90 = zeros(hor,n,numT);
SqHi68 = zeros(hor,n,numT);
SqLo68 = zeros(hor,n,numT);
SqHi90 = zeros(hor,n,numT);
SqLo90 = zeros(hor,n,numT);
TotHi68 = zeros(hor,n,numT);
TotLo68 = zeros(hor,n,numT);
TotHi90 = zeros(hor,n,numT);
TotLo90 = zeros(hor,n,numT);

% Fit the model separately for each GPR series
for ii = 1:numT

    % Pick the current GPR series and pair it with the macro data.
    idxT = find(strcmp(CountryNames, target_names{ii}));
    if isempty(idxT)
        error('Series %s not found in CountryNames.', target_names{ii});
    end
    GPR_series = NEW_28vars(:,idxT);
    vardata = [GPR_series Macro];

    % Extract the identified shock, then add its squared term to SVARX.
    [D,~,~,~,~,~,~,~,~,~,eta] = chol_irf(vardata,n,p,c,hor);
    finshock = eta(:,pos_shock);
    finshock_squared  = eta(:,pos_shock).^2;
    finaldata = vardata(p+1:end,:);
    news = [finshock finshock_squared];
    m = size(news,2);

    % SVARX supplies response estimates and confidence bounds.
    SVARX;

    % Scale responses to a two-standard-deviation shock.
    s = 2 * std(finshock);

    % Linear and quadratic terms occupy adjacent columns.
    for v_idx = 1:n
        col_lin = (v_idx-1)*m + 1;
        col_sq  = col_lin + 1;
        lin_irf = s    * C_wold(:, col_lin);
        lin_hi90  = s    * HighC90(:, col_lin);
        lin_lo90  = s    * LowC90(:, col_lin);
        lin_hi68  = s    * HighC(:,   col_lin);
        lin_lo68  = s    * LowC(:,    col_lin);
        sq_irf = (s^2)* C_wold(:, col_sq);
        sq_hi90 = (s^2)* HighC90(:, col_sq);
        sq_lo90 = (s^2)* LowC90(:, col_sq);
        sq_hi68 = (s^2)* HighC(:,   col_sq);
        sq_lo68 = (s^2)* LowC(:,    col_sq);

        % Combine both contributions and their confidence bounds.
        tot_irf = lin_irf + sq_irf;
        lin_up90 = lin_hi90 - lin_irf;
        lin_dn90 = lin_irf  - lin_lo90;
        sq_up90  = sq_hi90  - sq_irf;
        sq_dn90  = sq_irf   - sq_lo90;
        tot_up90 = sqrt(lin_up90.^2 + sq_up90.^2);
        tot_dn90 = sqrt(lin_dn90.^2 + sq_dn90.^2);
        tot_hi90 = tot_irf + tot_up90;
        tot_lo90 = tot_irf - tot_dn90;
        lin_up68 = lin_hi68 - lin_irf;
        lin_dn68 = lin_irf  - lin_lo68;
        sq_up68  = sq_hi68  - sq_irf;
        sq_dn68  = sq_irf   - sq_lo68;
        tot_up68 = sqrt(lin_up68.^2 + sq_up68.^2);
        tot_dn68 = sqrt(lin_dn68.^2 + sq_dn68.^2);
        tot_hi68 = tot_irf + tot_up68;
        tot_lo68 = tot_irf - tot_dn68;

        % Keep the three types of response for plotting.
        LinIRF(:,v_idx,ii)  = lin_irf;
        SqIRF(:,v_idx,ii) = sq_irf;
        TotIRF(:,v_idx,ii)  = tot_irf;
        LinHi68(:,v_idx,ii) = lin_hi68;
        LinLo68(:,v_idx,ii) = lin_lo68;
        LinHi90(:,v_idx,ii) = lin_hi90;
        LinLo90(:,v_idx,ii) = lin_lo90;
        SqHi68(:,v_idx,ii)  = sq_hi68;
        SqLo68(:,v_idx,ii)  = sq_lo68;
        SqHi90(:,v_idx,ii)  = sq_hi90;
        SqLo90(:,v_idx,ii)  = sq_lo90;
        TotHi68(:,v_idx,ii) = tot_hi68;
        TotLo68(:,v_idx,ii) = tot_lo68;
        TotHi90(:,v_idx,ii) = tot_hi90;
        TotLo90(:,v_idx,ii) = tot_lo90;
    end
end

% Linear response
plot_irf_3shocks( ...
    LinIRF, LinLo68, LinHi68, LinLo90, LinHi90, ...
    VARnames, tvec, target_names, target_colors, ...
    'Linear IRFs: GPR(World) vs GPR-Act vs GPR-Threat');

% Quadratic response
plot_irf_3shocks( ...
    SqIRF, SqLo68, SqHi68, SqLo90, SqHi90, ...
    VARnames, tvec, target_names, target_colors, ...
    'Quadratic IRFs: GPR(World) vs GPR-Act vs GPR-Threat');

% Total response
plot_irf_3shocks( ...
    TotIRF, TotLo68, TotHi68, TotLo90, TotHi90, ...
    VARnames, tvec, target_names, target_colors, ...
    'Total IRFs: GPR(World) vs GPR-Act vs GPR-Threat');

% Plot all nine responses for one IRF component.
function plot_irf_3shocks(IRF,Lo68,Hi68,Lo90,Hi90, ...
                           VARnames,tvec,target_names,target_colors,TitleName)
figure;
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
nvar = numel(VARnames);
numT = numel(target_names);

% Softer fills so overlapping intervals remain visible
CI_colors = [ ...
    0.50 0.50 0.50;
    1.00 0.40 0.40;
    0.30 0.40 1.00];
for v = 1:nvar
    nexttile; hold on; box on;
    plot(tvec,0*tvec,'k-','LineWidth',0.5,'HandleVisibility','off');
    for k = 1:numT
        fill([tvec fliplr(tvec)], ...
             [Lo90(:,v,k)' fliplr(Hi90(:,v,k)')], ...
             CI_colors(k,:), ...
             'FaceAlpha',0.12, 'EdgeColor','none', ...
             'HandleVisibility','off');
        fill([tvec fliplr(tvec)], ...
             [Lo68(:,v,k)' fliplr(Hi68(:,v,k)')], ...
             CI_colors(k,:), ...
             'FaceAlpha',0.22, 'EdgeColor','none', ...
             'HandleVisibility','off');
        plot(tvec, IRF(:,v,k), ...
             'Color', target_colors(k,:), 'LineWidth',2.5);
    end
    title(VARnames{v}, 'FontSize', 16, 'FontWeight', 'bold');
    grid on; xlim([tvec(1) tvec(end)]);
    set(gca,'FontSize',12);
    if v == 1
        legend(target_names,'NumColumns',3,'Location','southoutside',...
               'FontSize',12);
    end
end
sgtitle(TitleName,'FontSize',18,'FontWeight','bold');
end
