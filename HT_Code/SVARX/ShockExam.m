% Cross-country robustness checks, quarterly sample (1994Q2–2023Q2)
% Cholesky-identified GPR shocks; no LASSO step.

clear; clc;
addpath('Functions');

% Load the aligned quarterly data
load('AlignedData_1994Q2_2023Q2.mat','NEW_28vars','NNEW');
NEW_9 = NNEW;
NEW_28Q = NEW_28vars;
[T_9, K_9] = size(NEW_9);
[T_28, K_28] = size(NEW_28Q);
if T_9 ~= T_28
    error('Length mismatch: NNEW and NEW_28vars must have same T.');
end
T = T_9;
numSeries = K_28;

% Order of the 28 GPR series
CountryNames = { ...
    'GPR(World)','GPR-Threat','GPR-Act', ...
    'Australia','Brazil','Canada','Switzerland','China','Germany','Egypt', ...
    'France','United Kingdom','India','Israel','Japan','South Korea', ...
    'Mexico','Netherlands','Russia','Saudi Arabia','Sweden','Turkey', ...
    'Taiwan','Ukraine','United States','Venezuela','Vietnam','South Africa' ...
};
if numel(CountryNames) ~= numSeries
    error('CountryNames length (%d) ≠ numSeries (%d).',numel(CountryNames),numSeries);
end

% VAR and test settings
p = 1;
c = 1;
hor = 1;
pos_shock = 1;
nwLag = 1;
maxLagLB  = 12;
archLag = 4;
Results = cell(numSeries,1);
LB_Results = struct();
Norm_Results = struct();
ARCH_Results = struct();

for j = 1:numSeries
    fprintf('\n-------------------------------------------------\n');
    fprintf(' Processing %2d / %2d : %s\n', j, numSeries, CountryNames{j});
    fprintf('-------------------------------------------------\n');

    % Swap the GPR column; the eight macro series stay the same.
    tmp = NEW_9;
    tmp(:,1) = NEW_28Q(:,j);
    GPR = tmp(:,1);
    controls = tmp(:,2:end);
    vardata  = tmp(:,1:9);
    R = struct();
    R.Country = CountryNames{j};

    % Check whether the fitted VAR is stable.
    R.Stability = struct('Eigenvalues',[],'IsStable',NaN);
    try
        [Acomp, SIGMA_var, U_var, V_var] = var_stability(vardata, p);
        eigvals  = eig(Acomp);
        isStable = all(abs(eigvals) < 1);
        R.Stability.Eigenvalues = eigvals;
        R.Stability.IsStable = isStable;
        if isStable
            fprintf('  VAR Stable ✓ (all eigenvalues < 1)\n');
        else
            fprintf('  VAR Unstable ✗ (some eigenvalues ≥ 1)\n');
        end
    catch ME
        fprintf('  Stability test FAILED: %s\n', ME.message);
    end

    % The first Cholesky shock is the GPR shock.
    [Tvar, nvar] = size(vardata);
    n = nvar;
    try
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, eta] = chol_irf(vardata, n, p, c, hor);
    catch ME
        error('chol_irf failed for %s: %s', CountryNames{j}, ME.message);
    end
    shock = eta(:, pos_shock);
    T_eta = length(shock);
    R.Eta = shock;

    % Match the shock at t with the next quarter and current controls.
    eta_t = shock(1:end-1);
    eta_tp1 = shock(2:end);
    GPR_t_reg = vardata(p+1:end-1, 1);
    X_t_reg = vardata(p+1:end-1, 2:end);
    Treg = length(eta_tp1);
    if any([length(eta_t), length(GPR_t_reg), size(X_t_reg,1)] ~= Treg)
        error('Length mismatch in regression vars for j = %d.', j);
    end
    Xreg = [ones(Treg,1), eta_t, GPR_t_reg, X_t_reg];

    % Predictability regression with Newey-West standard errors.
    [beta, se, tstat, pval] = newey_west(Xreg, eta_tp1, nwLag);
    R.Beta  = beta;
    R.SE = se;
    R.P = pval;
    R.Tstat = tstat;
    Results{j} = R;
end

% Serial correlation in the extracted shocks
for j = 1:numSeries
    eta_j = Results{j}.Eta;
    [h,pval,stat] = lbqtest(eta_j, 'Lags', maxLagLB);
    LB_Results(j).Country = Results{j}.Country;
    LB_Results(j).LB_stat = stat;
    LB_Results(j).LB_pval = pval;
    LB_Results(j).Reject  = h;
    fprintf('%-15s  LBQ stat = %8.3f   p = %6.4f   reject? %d\n', ...
        LB_Results(j).Country, stat, pval, h);
end

% Normality and distribution shape
for j = 1:numSeries
    eta_j = Results{j}.Eta;
    [h_jb, p_jb, jbstat] = jbtest(eta_j);
    sk = skewness(eta_j);
    kt = kurtosis(eta_j);
    Norm_Results(j).Country  = Results{j}.Country;
    Norm_Results(j).JB_stat  = jbstat;
    Norm_Results(j).JB_pval  = p_jb;
    Norm_Results(j).RejectJB = h_jb;
    Norm_Results(j).Skewness = sk;
    Norm_Results(j).Kurtosis = kt;
    fprintf('%-15s  JB = %8.3f   p = %6.4f   skew = %6.3f   kurt = %6.3f   reject? %d\n', ...
        Norm_Results(j).Country, jbstat, p_jb, sk, kt, h_jb);
end

% ARCH effects in the shocks
for j = 1:numSeries
    eta_j = Results{j}.Eta;
    [h_arch, p_arch, stat_arch] = archtest(eta_j, 'Lags', archLag);
    ARCH_Results(j).Country = Results{j}.Country;
    ARCH_Results(j).ARCH_stat = stat_arch;
    ARCH_Results(j).ARCH_pval = p_arch;
    ARCH_Results(j).Reject = h_arch;
    fprintf('%-15s  ARCH LM = %8.3f   p = %6.4f   reject? %d\n', ...
        ARCH_Results(j).Country, stat_arch, p_arch, h_arch);
end

% VAR stability table
fprintf('\n\n%% ==========================================\n');
fprintf('%% LaTeX Table: VAR Stability Across Countries\n');
fprintf('%% ==========================================\n');
fprintf('\\begin{table}[h!]\n');
fprintf('\\centering\n');
fprintf('\\caption{VAR Stability Test Across Countries (Eigenvalues $< 1$)}\n');
fprintf('\\begin{tabular}{lcc}\n');
fprintf('\\hline\\hline\n');
fprintf('Country & Stable? & Max $|\\lambda|$ \\\\\\n');
fprintf('\\hline\n');

for j = 1:numSeries
    eigvals = Results{j}.Stability.Eigenvalues;
    if isempty(eigvals)
        fprintf('%s & NA & NA \\\\\\n', Results{j}.Country);
    else
        maxeig = max(abs(eigvals));
        stab = Results{j}.Stability.IsStable;
        fprintf('%s & %d & %6.3f \\\\\\n', ...
            Results{j}.Country, stab, maxeig);
    end
end
fprintf('\\hline\\hline\n');
fprintf('\\end{tabular}\n');
fprintf('\\end{table}\n');

% Full predictability regression table
fprintf('\n\n%% ==========================================\n');
fprintf('%% LaTeX Table: Full Predictability Regression (Cholesky)\n');
fprintf('%% ==========================================\n');
fprintf('\\begin{sidewaystable}[p]\n');
fprintf('\\scriptsize\n');
fprintf('\\centering\n');
fprintf(['\\caption{Predictability Regression Coefficients for Cholesky-Identified ' ...
         'GPR Structural Shocks Across Countries} \n']);
fprintf('\\label{tab:chol_predictability_full}\n');
fprintf(['\\begin{tabular}{l' repmat('cc',1,10) '}\n']);
fprintf('\\hline\\hline\n');
fprintf(['Country & $\\beta_1$ & p($\\beta_1$) & $\\beta_2$ & p($\\beta_2$) ' ...
         '& $\\gamma_1$ & p($\\gamma_1$) & $\\gamma_2$ & p($\\gamma_2$) ' ...
         '& $\\gamma_3$ & p($\\gamma_3$) & $\\gamma_4$ & p($\\gamma_4$) ' ...
         '& $\\gamma_5$ & p($\\gamma_5$) & $\\gamma_6$ & p($\\gamma_6$) ' ...
         '& $\\gamma_7$ & p($\\gamma_7$) & $\\gamma_8$ & p($\\gamma_8$) \\\\\\n']);
fprintf('\\hline\n');

for j = 1:numSeries
    b = Results{j}.Beta;
    p = Results{j}.P;
    beta1  = b(2);  p1 = p(2);
    beta2  = b(3);  p2 = p(3);
    gamma  = b(4:end);
    pgamma = p(4:end);
    fprintf(['%s & %7.3f & %7.3f & %7.3f & %7.3f ' ...
             '& %7.3f & %7.3f & %7.3f & %7.3f ' ...
             '& %7.3f & %7.3f & %7.3f & %7.3f ' ...
             '& %7.3f & %7.3f & %7.3f & %7.3f ' ...
             '& %7.3f & %7.3f & %7.3f & %7.3f \\\\\\n'], ...
             Results{j}.Country, ...
             beta1, p1, beta2, p2, ...
             gamma(1), pgamma(1), gamma(2), pgamma(2), ...
             gamma(3), pgamma(3), gamma(4), pgamma(4), ...
             gamma(5), pgamma(5), gamma(6), pgamma(6), ...
             gamma(7), pgamma(7), gamma(8), pgamma(8));
end
fprintf('\\hline\\hline\n');
fprintf('\\end{tabular}\n');
fprintf('\\end{sidewaystable}\n');

% Ljung-Box table
fprintf('\n\n%% ==========================================\n');
fprintf('%% LaTeX Table: Ljung--Box Q-test for $\\eta_t$ Shocks (Quarterly)\n');
fprintf('%% ==========================================\n');
fprintf('\\begin{table}[h!]\\centering\n');
fprintf('\\begin{tabular}{lcc}\\hline\\hline\n');
fprintf('Country & Q-statistic & p-value \\\\\\n');
fprintf('\\hline\n');

for j = 1:numSeries
    fprintf('%s & %7.3f & %7.4f \\\\\\n', ...
        LB_Results(j).Country, LB_Results(j).LB_stat, LB_Results(j).LB_pval);
end
fprintf('\\hline\\hline\\end{tabular}\\end{table}\n');

% Normality table
fprintf('\n\n%% ==========================================\n');
fprintf('%% LaTeX Table: Normality (Jarque--Bera) for $\\eta_t$ Shocks\n');
fprintf('%% ==========================================\n');
fprintf('\\begin{table}[h!]\\centering\n');
fprintf('\\begin{tabular}{lcccc}\\hline\\hline\n');
fprintf('Country & JB-statistic & p-value & Skewness & Kurtosis \\\\\\n');
fprintf('\\hline\n');

for j = 1:numSeries
    fprintf('%s & %7.3f & %7.4f & %7.3f & %7.3f \\\\\\n', ...
        Norm_Results(j).Country, ...
        Norm_Results(j).JB_stat, ...
        Norm_Results(j).JB_pval, ...
        Norm_Results(j).Skewness, ...
        Norm_Results(j).Kurtosis);
end
fprintf('\\hline\\hline\\end{tabular}\\end{table}\n');

% ARCH table
fprintf('\n\n%% ==========================================\n');
fprintf('%% LaTeX Table: ARCH LM Test for $\\eta_t$ Shocks\n');
fprintf('%% ==========================================\n');
fprintf('\\begin{table}[h!]\\centering\n');
fprintf('\\begin{tabular}{lcc}\\hline\\hline\n');
fprintf('Country & ARCH-statistic & p-value \\\\\\n');
fprintf('\\hline\n');

for j = 1:numSeries
    fprintf('%s & %7.3f & %7.4f \\\\\\n', ...
        ARCH_Results(j).Country, ...
        ARCH_Results(j).ARCH_stat, ...
        ARCH_Results(j).ARCH_pval);
end
fprintf('\\hline\\hline\\end{tabular}\\end{table}\n');

% Summary table
fprintf('\n\n%% ==========================================\n');
fprintf('%% Combined Robustness Table (Quarterly, Cholesky shocks)\n');
fprintf('%% ==========================================\n\n');
fprintf('\\begin{table}[H]\n');
fprintf('\\centering\n');
fprintf(['\\caption{Robustness Checks for Predictability of ' ...
         'Cholesky-Identified GPR Structural Shocks (Quarterly Sample 1994Q2--2023Q2)}\n']);
fprintf('\\label{tab:gpr_robustness_combined_Q}\n');
fprintf(['\\begin{tabular}{lcccccc}\n' ...
         '\\hline\\hline\n' ...
         'Country & $\\beta_2$ & p-value($\\beta_2$) & Q(12) p-value & JB p-value & ARCH p-value & Stable VAR? \\\\\\n' ...
         '\\hline\n']);

for j = 1:numSeries
    b2 = Results{j}.Beta(3);
    p2 = Results{j}.P(3);
    pLB = LB_Results(j).LB_pval;
    pJB = Norm_Results(j).JB_pval;
    pARCH = ARCH_Results(j).ARCH_pval;
    stab = Results{j}.Stability.IsStable;
    fprintf('%s & %7.3f & %7.3f & %7.4f & %7.4f & %7.4f & %d \\\\\\n', ...
        Results{j}.Country, b2, p2, pLB, pJB, pARCH, stab);
end
fprintf('\\hline\\hline\n');
fprintf('\\end{tabular}\n');
fprintf(['\\begin{flushleft}\\footnotesize ', ...
         'Notes: Columns report (i) Newey--West estimate of $\\beta_2$ (coefficient on GPR$_t$) and its p-value, ', ...
         '(ii) Ljung--Box Q(12) p-value for the Cholesky-identified shock $\\eta_t$, ', ...
         '(iii) Jarque--Bera normality test p-value for $\\eta_t$, ', ...
         '(iv) ARCH LM test p-value for $\\eta_t$ (up to lag ', num2str(archLag), '), ', ...
         'and (v) VAR stability indicator (1 = all eigenvalues inside unit circle). ', ...
         'Sample: 1994Q2--2023Q2. \\end{flushleft}\n']);
fprintf('\\end{table}\n');
