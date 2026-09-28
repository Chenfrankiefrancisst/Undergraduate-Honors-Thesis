clear; clc; close all;

addpath("ADLdata")

% Start by a random country: say, US
NEWUS2;   

%% Parameters
input   = 2;
t0      = 1;               
A       = 16;   % forecast horizon
L_imp   = A;
level   = 0;
alpha   = 0.68;
N_B     = 10000;
mindelay= 0;
IC      = 0;
trend   = 0;
h       = 0:L_imp;
lw      = 2;
fs      = 20;
scale = exist('input','var') * input + ~exist('input','var'); 
xx      = [0:1:L_imp, L_imp:-1:0];

% variables of interest
colNames = {
    'CPI', ...
    'US Policy Rate', ...
    'One Year Yield', ...
    'Gold Demand', ...
    'Country-Specific Gold Production', ...
    'Real Economic Activity', ...
    'Gold Spot Price', ...
    'Durable' ...
};

nVars = size(NEW, 2) - 1;  

% iterate among countries:
dataScripts   = {'NEWW2', 'NEWGPA2', 'NEWGPT2', 'NEWUS2', 'NEWCA2', 'NEWAU2', 'NEWME2', 'NEWCN2', 'NEWRU2'};
countryLabels = {'W',     'WThreat',       'WAct',       'US',     'CA',     'AU',     'ME',     'CN',     'RU'};

%%  main iteration

for j = 2:(nVars+1)   
    
    if (j-1) <= numel(colNames)
        varName = colNames{j-1};
    else
        varName = sprintf('Var %d', j-1);
    end
    
    figName = sprintf('IRFs of %s to GPR shocks', varName);
    figure('Name', figName, 'Position', [50, 50, 1200, 800]);
    
    for c = 1:numel(countryLabels)
        
        eval(dataScripts{c});   
        
        X  = diff(NEW(:, 1));
        %X  = NEW(:, 1) ;
        t1 = length(X);         
        
        y = NEW(t0:t1, j);
        
        [impBIC, cbBIC, ~, lag_length, ~] = ir_ADL( ...
            y, X, 0:10, A, level, L_imp, alpha, N_B, mindelay, IC, trend);
        
        imp = impBIC * scale;
        if size(cbBIC,1) == 2
            lo = cbBIC(1,:) * scale;
            hi = cbBIC(2,:) * scale;
        else
            lo = cbBIC(:,1)' * scale;
            hi = cbBIC(:,2)' * scale;
        end
        
        yy = [lo, fliplr(hi)];
        
        subplot(3, 3, c);
        hold on;
        
        fill(xx, yy, [0.8 0.8 0.8], 'EdgeColor','none');   
        plot(h, lo, 'k--', 'LineWidth', 1);
        plot(h, hi, 'k--', 'LineWidth', 1);
        plot(h, imp, 'k-',  'LineWidth', lw);
        yline(0, 'r:'); 
        
        xlabel('Quarter', 'Interpreter', 'latex', 'FontSize', fs-4); 
        ylabel('$\psi_{x,c}$', 'Interpreter', 'latex', 'FontSize', fs-4);
        
        thisTitle = sprintf('%s - %s (B$_{opt}$ = %d)', ...
                            countryLabels{c}, varName, lag_length.I);
        title(thisTitle, 'Interpreter', 'latex', 'FontSize', fs-3);
        
        xlim([h(1) h(end)]);
        set(gca, 'Layer', 'top', 'TickLabelInterpreter', 'latex', 'FontSize', fs-6);
        
        if c == 1
            legend({'68\% band','band edge','band edge','IRF','zero'}, ...
                'Location','best', 'FontSize', 10, 'Interpreter', 'latex');
        end
        
        hold off;
    end
    
    sgtitle(sprintf('Impulse Responses of %s to GPR shocks (US, CA, AU, ME, CN, RU)', varName), ...
        'FontSize', fs, 'Interpreter', 'latex');
end
