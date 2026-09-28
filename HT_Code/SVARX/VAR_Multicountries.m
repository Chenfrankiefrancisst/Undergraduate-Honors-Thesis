% Cross-country GPR responses: linear, quadratic, and total IRFs
clear; clc; close all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
rng(0);
addpath('Functions');
addpath('Data');

% Load the aligned GPR and macro series
load('AlignedData_1994Q2_2023Q2.mat','NEW_28vars','NNEW')
Macro = NNEW(:,2:9);
T = size(Macro,1);
n = 9;

% Country order matches the columns of NEW_28vars
CountryNames = { ...
    'GPR(World)','GPR-Threat','GPR-Act', ...
    'Australia','Brazil','Canada','Switzerland','China','Germany','Egypt', ...
    'France','United Kingdom','India','Israel','Japan','South Korea', ...
    'Mexico','Netherlands','Russia','Saudi Arabia','Sweden','Turkey', ...
    'Taiwan','Ukraine','United States','Venezuela','Vietnam','South Africa' ...
};

% Labels used in the nine panels
VARnames = { ...
    'GPR Variable';
    'CPIU';
    'Shadow Rate';
    '1Y Gov Yield';
    'Gold Demand';
    'Gold Production';
    'Real Activity';
    'Gold Spot';
    'Durables';
};

% SVAR settings
p   = 1;
c   = 1;
hor = 49;
cum = 0;
pos_shock = 1;
nonlin = 1;
restr  = 2;
m      = 2;
q      = 0;
tvec = 0:hor-1;

% Country comparisons to plot
Groups = {
 {'GPR(World)','GPR-Threat','GPR-Act','Russia','Ukraine','China'};
 {'GPR(World)','GPR-Threat', 'GPR-Act','United States','United Kingdom','France','Canada'};
 {'GPR(World)','GPR-Threat','China','Japan','South Korea'};
 {'GPR(World)','GPR-Threat','China','Taiwan','South Korea'};
 {'GPR(World)','GPR-Threat','GPR-Act','Russia','Ukraine'};
 {'GPR(World)','United States','Japan','South Korea'};
 {'GPR(World)','GPR-Act','GPR-Threat','United States','Israel'};
 {'GPR(World)','United States','Saudi Arabia','Turkey','Egypt','Israel'};
};
numGroups = numel(Groups);
CI_color = [160 160 160]/255;

% Estimate and plot each group
for G = 1:numGroups
    current_group = Groups{G};
    numCountries  = numel(current_group);
    cmap = dynamic_colors(numCountries);
    idx = zeros(numCountries,1);
    for i=1:numCountries
        idx(i) = find(strcmp(CountryNames,current_group{i}));
    end
    LinIRF = zeros(hor,n,numCountries);
    SqIRF  = zeros(hor,n,numCountries);
    TotIRF = zeros(hor,n,numCountries);
    TotHigh90 = zeros(hor,n);
    TotLow90  = zeros(hor,n);

    % Estimate the world series first; its intervals are used in the plots.
    world_idx = find(strcmp(current_group,'GPR(World)'));
    GPR_series = NEW_28vars(:, idx(world_idx));

    % Pair the current GPR measure with the same macro series.
    vardata = [GPR_series Macro];
    [~,~,~,~,~,~,~,~,~,~,eta] = chol_irf(vardata,n,p,c,hor);
    fin = eta(:,1);
    fin2 = fin.^2;
    finaldata = vardata(p+1:end,:);
    news      = [fin fin2];

    % SVARX writes the response estimates and confidence bounds.
    SVARX;
    s = 2 * std(fin);
    for v=1:n
        col_lin = (v-1)*2 + 1;
        col_sq  = col_lin + 1;
        lin = s    * C_wold(:,col_lin);
        sq  = (s^2)* C_wold(:,col_sq);
        tot = lin + sq;
        LinIRF(:,v,world_idx) = lin;
        SqIRF(:,v,world_idx)  = sq;
        TotIRF(:,v,world_idx) = tot;
        lin_hi = s    * HighC90(:,col_lin);
        lin_lo = s    * LowC90(:,col_lin);
        sq_hi  = (s^2)* HighC90(:,col_sq);
        sq_lo  = (s^2)* LowC90(:,col_sq);
        up = sqrt((lin_hi-lin).^2 + (sq_hi-sq).^2);
        dn = sqrt((lin-lin_lo).^2 + (sq-sq_lo).^2);

        TotHigh90(:,v) = tot + up;
        TotLow90(:,v)  = tot - dn;
    end

    % Run the same calculation for the remaining countries.
    for ctry = 1:numCountries
        if ctry == world_idx, continue; end
        GPR_series = NEW_28vars(:, idx(ctry));

        % Use the same macro series for this country.
        vardata = [GPR_series Macro];
        [~,~,~,~,~,~,~,~,~,~,eta] = chol_irf(vardata,n,p,c,hor);
        fin = eta(:,1);
        fin2 = fin.^2;
        finaldata = vardata(p+1:end,:);
        news      = [fin fin2];

        % Update the responses for this shock.
        SVARX;
        for v=1:n
            col_lin = (v-1)*2 + 1;
            col_sq  = col_lin + 1;
            lin = s    * C_wold(:,col_lin);
            sq  = (s^2)* C_wold(:,col_sq);
            tot = lin + sq;
            LinIRF(:,v,ctry) = lin;
            SqIRF(:,v,ctry)  = sq;
            TotIRF(:,v,ctry) = tot;
        end
    end

% Linear responses
figure;
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
for v=1:n
    nexttile; hold on; box on;
    plot(tvec,0*tvec,'k-','LineWidth',0.3,'HandleVisibility','off');
    hi = TotHigh90(:,v) - SqIRF(:,v,world_idx);
    lo = TotLow90(:,v)  - SqIRF(:,v,world_idx);
    fill([tvec fliplr(tvec)], [lo' fliplr(hi')], CI_color, ...
        'FaceAlpha',.22,'EdgeColor','none','HandleVisibility','off');
    h = gobjects(numCountries,1);
    h(world_idx) = plot(tvec,LinIRF(:,v,world_idx), ...
        'Color',cmap(world_idx,:),'LineWidth',2.5);

    % Run the same calculation for the remaining countries.
    for ctry = 1:numCountries
        if ctry==world_idx, continue; end
        h(ctry) = plot(tvec,LinIRF(:,v,ctry), ...
            'Color',cmap(ctry,:),'LineWidth',2);
    end
    title(VARnames{v});
    grid on; set(gca,'FontSize',10); xlim([0 hor-1]);
    if v==1
        L = current_group;
        L{world_idx} = [L{world_idx} ' (CI)'];
        legend(h,L,'NumColumns',numCountries);
    end
end
sgtitle(['Linear IRFs - Group ' num2str(G)]);

% Quadratic responses
figure;
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
for v=1:n
    nexttile; hold on; box on;
    plot(tvec,0*tvec,'k-','LineWidth',0.3,'HandleVisibility','off');
    hi = TotHigh90(:,v) - LinIRF(:,v,world_idx);
    lo = TotLow90(:,v)  - LinIRF(:,v,world_idx);
    fill([tvec fliplr(tvec)], [lo' fliplr(hi')], CI_color, ...
        'FaceAlpha',.22,'EdgeColor','none','HandleVisibility','off');
    h = gobjects(numCountries,1);
    h(world_idx) = plot(tvec,SqIRF(:,v,world_idx), ...
        'Color',cmap(world_idx,:),'LineWidth',2.5);

    % Run the same calculation for the remaining countries.
    for ctry = 1:numCountries
        if ctry==world_idx, continue; end
        h(ctry) = plot(tvec,SqIRF(:,v,ctry), ...
            'Color',cmap(ctry,:),'LineWidth',2);
    end
    title(VARnames{v});
    grid on; xlim([0 hor-1]);
    if v==1
        L = current_group;
        L{world_idx} = [L{world_idx} ' (CI)'];
        legend(h,L,'NumColumns',numCountries);
    end
end
sgtitle(['Quadratic IRFs - Group ' num2str(G)]);

% Total responses
figure;
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
for v=1:n
    nexttile; hold on; box on;
    plot(tvec,0*tvec,'k-','LineWidth',0.3,'HandleVisibility','off');
    hi = TotHigh90(:,v);
    lo = TotLow90(:,v);
    fill([tvec fliplr(tvec)], [lo' fliplr(hi')], CI_color, ...
        'FaceAlpha',.22,'EdgeColor','none','HandleVisibility','off');
    h = gobjects(numCountries,1);
    h(world_idx) = plot(tvec,TotIRF(:,v,world_idx), ...
        'Color',cmap(world_idx,:),'LineWidth',2.5);

    % Run the same calculation for the remaining countries.
    for ctry = 1:numCountries
        if ctry==world_idx, continue; end
        h(ctry) = plot(tvec,TotIRF(:,v,ctry), ...
            'Color',cmap(ctry,:),'LineWidth',2);
    end
    title(VARnames{v});
    grid on; xlim([0 hor-1]);
    if v==1
        L = current_group;
        L{world_idx} = [L{world_idx} ' (CI)'];
        legend(h,L,'NumColumns',numCountries);
    end
end
sgtitle(['Total IRFs - Group ' num2str(G)]);
end

% Colors for the country lines
function C = dynamic_colors(numCountries)
    base = [
        0   0   0;
        180 30 30;
        20 70 200;
        0 120 0;
        140 20 160;
        200 100 0;
        120 70 40;
        200 0 120;
        40 150 200;
        0 160 160;
        100 0 0;
        0 0 120;
    ] / 255;
    if numCountries <= size(base,1)
        C = base(1:numCountries,:);
    else
        k = numCountries - size(base,1);
        hsv = [linspace(0,1,k)' ones(k,1)*0.8 ones(k,1)*0.9];
        extra = hsv2rgb(hsv);
        C = [base; extra];
    end
end
