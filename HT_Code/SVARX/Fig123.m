clear; clc; close all;

%% Data loading
GPR_All; 

% original order: just for reference
% 1  GPR
% 2  GPRT
% 3  GPRA
% 4  Australia
% 5  Brazil
% 6  Canada
% 7  Switzerland
% 8  China
% 9  Germany
% 10 Egypt
% 11 France
% 12 United Kingdom
% 13 India
% 14 Israel
% 15 Japan
% 16 South Korea
% 17 Mexico
% 18 Netherlands
% 19 Russia
% 20 Saudi Arabia
% 21 Sweden
% 22 Turkey
% 23 Taiwan
% 24 Ukraine
% 25 United States
% 26 Venezuela
% 27 Vietnam
% 28 South Africa

%% reordering

newOrder = [ ...
    1, 2, 3, ...                           % GPR, GPRT, GPRA
    25, 6, 4, 12, 11, 9, 18, 21, 7, ...    % US, Canada, Australia, UK, France, Germany, Netherlands, Sweden, Switzerland
    15, 16, 8, 23, 27, 13, 17, 5, ...      % Japan, Korea, China, Taiwan, Vietnam, India, Mexico, Brazil
    19, 24, 14, 20, 22, 10, 26, 28];       % Russia, Ukraine, Israel, Saudi, Turkey, Egypt, Venezuela, South Africa

NEW_reordered = NEW(:, newOrder);

newNames = { ...
    'GPR', 'GPRT', 'GPRA', ...
    'United States', 'Canada', 'Australia', ...
    'United Kingdom', 'France', 'Germany', 'Netherlands', 'Sweden', 'Switzerland', ...
    'Japan', 'South Korea', 'China', 'Taiwan', 'Vietnam', 'India', 'Mexico', 'Brazil', ...
    'Russia', 'Ukraine', 'Israel', 'Saudi Arabia', 'Turkey', 'Egypt', 'Venezuela', 'South Africa'};

%%  correlation matrix

CorrMatrix = corr(NEW_reordered, 'rows', 'pairwise');

%% heatmap

figure('Name','Correlation Matrix of GPR Across Countries','Position',[100 100 1000 800]);

h = heatmap(newNames, newNames, CorrMatrix, ...
    'Colormap', turbo, ...          
    'ColorLimits', [-1 1], ...
    'CellLabelFormat','%.2f');

h.FontName = 'Palatino';
h.FontSize = 14;
h.XLabel = 'GPR / Country / Region';
h.YLabel = 'GPR / Country / Region';
h.Title = 'Correlation of GPR Measures and Country GPR Indices';

h.XDisplayLabels = strrep(h.XDisplayLabels, ' ', '\_'); 



%% Take 63 quarters（2010Q1–2025Q3) for recent geopolitical changes
T = 63;   

GPRmat_last = NEW_reordered(end-T+1:end, :);   

CorrMatrix = corr(GPRmat_last, 'rows','pairwise');

figure('Name','Correlation Matrix of GPR (Last 63 Quarters)', ...
       'Position',[100 100 1100 850]);

h = heatmap(newNames, newNames, CorrMatrix, ...
    'Colormap', turbo, ...         
    'ColorLimits', [-1 1], ...
    'CellLabelFormat','%.2f');

h.FontName = 'Palatino';
h.FontSize = 14;

h.Title = sprintf('Correlation of GPR (Last %d Quarters: 2010Q1–2025Q3)', T);



%%  Define names in EXACT column order 
countryNames_exact = { ...
    'GPR', 'GPRT', 'GPRA', ...
    'Australia', 'Brazil', 'Canada', 'Switzerland', 'China', 'Germany', ...
    'Egypt', 'France', 'United Kingdom', 'India', 'Israel', ...
    'Japan', 'South Korea', 'Mexico', 'Netherlands', ...
    'Russia', 'Saudi Arabia', 'Sweden', 'Turkey', ...
    'Taiwan', 'Ukraine', 'United States', ...
    'Venezuela', 'Vietnam', 'South Africa'};

std_all = std(NEW, 0, 1, 'omitnan');    % 1 × 28 vector

[std_sorted, idx] = sort(std_all);
names_sorted = countryNames_exact(idx);

figure('Position',[300 300 1500 600]);

bar(std_sorted, 'FaceColor', [0.75 0.75 0.75]); % grey

set(gca, 'FontName','Palatino', ...
         'TickLabelInterpreter','latex', ...
         'FontSize',14);

xticks(1:length(std_sorted));
xticklabels(names_sorted);
xtickangle(45);

ylabel('Standard Deviation', ...
       'Interpreter','latex', ...
       'FontName','Palatino', ...
       'FontSize',16);

title('Standard Deviation of GPR Variables (Sorted)', ...
       'Interpreter','latex', ...
       'FontName','Palatino', ...
       'FontSize',18);

grid on

for i = 1:length(std_sorted)
    text(i, std_sorted(i) + 0.02 * max(std_sorted), ...
        sprintf('$%.3f$', std_sorted(i)), ...
        'Interpreter','latex', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',12, ...
        'FontName','Palatino');
end



%% Case 1: Full Sample 
Corr_full = corr(NEW_reordered, 'rows','pairwise');

disp('===================================================');
disp('   FULL SAMPLE CORRELATION MATRIX (ALL OBS)');
disp('===================================================');

CorrTable_full = array2table(Corr_full, ...
    'VariableNames', newNames, ...
    'RowNames', newNames);

disp(CorrTable_full);


%% Case 2: Last 63 Quarters
T = 63;
GPR_last = NEW_reordered(end-T+1:end, :);
Corr_last = corr(GPR_last, 'rows','pairwise');

disp(' ');
disp('===================================================');
disp(sprintf('   CORRELATION MATRIX (Last %d Quarters: 2010Q1–2025Q3)', T));
disp('===================================================');

CorrTable_last = array2table(Corr_last, ...
    'VariableNames', newNames, ...
    'RowNames', newNames);

disp(CorrTable_last);
