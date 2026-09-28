

function irf_plot_VARX_separate(MiddleD,HighD,LowD,HighD68,LowD68, ...
                                 VARnames,Shocknames_merged, ...
                                 colorBNDS,colorBNDS68, ...
                                 hor,n)
% irf_plot_VARX_threegrid
%
% Makes THREE separate figures.
% Each figure is a 3x3 grid (n=9 variables).
%
% Figure 1: Shocknames_merged{1}  (e.g. 'GPR Shock')
% Figure 2: Shocknames_merged{2}  (e.g. 'GPR Shock Squared')
% Figure 3: Shocknames_merged{3}  (e.g. 'Overall (pos vs neg)')
%
% Inputs:
%   MiddleD      [hor x (n*3)]      main IRF per shock block (black solid line)
%   HighD,LowD   [hor x (n*3)]      90% band envelopes
%   HighD68,LowD68 [hor x (n*3)]    68% band envelopes
%   VARnames     {n x 1}            variable names
%   Shocknames_merged {3 x 1}       titles for each shock figure
%   colorBNDS    [1x3]              light gray for 90% band
%   colorBNDS68  [1x3]              darker gray for 68% band
%   hor          scalar             # horizons
%   n            scalar             # variables (should be 9)

    set(0,'defaultAxesFontName','Times');
    set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1.5)

    num_shocks_to_plot = 3;  % column groups per variable

    for shock_idx = 1:num_shocks_to_plot

        figure('Color','w');
        sgtitle(Shocknames_merged{shock_idx}, 'FontSize', 18, 'Interpreter','Latex');

        for var_idx = 1:n

            subplot(3,3,var_idx);

            % column index for this (variable, this shock)
            % layout: [var1_sh1 var1_sh2 var1_sh3 | var2_sh1 ...]
            col = (var_idx-1)*num_shocks_to_plot + shock_idx;

            % --- 90% band (light gray fill)
            s1 = fill([0:hor-1, fliplr(0:hor-1)]', ...
                      [ HighD(:,col); flipud(LowD(:,col)) ], ...
                      colorBNDS, 'EdgeColor','none'); hold on;

            % --- 68% band (darker gray fill)
            s2 = fill([0:hor-1, fliplr(0:hor-1)]', ...
                      [ HighD68(:,col); flipud(LowD68(:,col)) ], ...
                      colorBNDS68, 'EdgeColor','none'); hold on;

            % --- IRF line (black solid)
            p1 = plot(0:hor-1, MiddleD(:,col), ...
                      'k','LineWidth',2.0); hold on;

            % zero line
            plot(get(gca,'XLim'), [0 0], 'k--','LineWidth',1.0);

            % panel title = variable name
            var_label = VARnames{var_idx};
            var_label = strrep(var_label,'_','\_');
            title(var_label,'FontSize',11,'Interpreter','Latex');

            xlim([0 hor-1]);
            axis tight;
            set(gca,'FontSize',10);

            % x-label only for bottom row (subplots 7,8,9)
            if var_idx >= 7
                xlabel('Horizon','FontSize',10)
            end
        end

        % one legend per figure
        han = axes('Position',[0 0 1 1],'Visible','off');
        legend(han,[s1,s2,p1], ...
            {'90% bands','68% bands','IRF'}, ...
            'Orientation','Horizontal', ...
            'FontSize',12, ...
            'Interpreter','none', ...
            'Location','northoutside');

    end
end
