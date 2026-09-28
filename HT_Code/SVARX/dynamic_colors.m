%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DYNAMIC COLOR MAP (UNLIMITED), WORLD ALWAYS BLACK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = dynamic_colors(numCountries)
    % World always black
    base = [
        0   0   0;         % 1. World (fixed black)
        180 30 30;         % 2. red
        20 70 200;         % 3. blue
        0 120 0;           % 4. green
        140 20 160;        % 5. purple
        200 100 0;         % 6. orange
        120 70 40;         % 7. brown
        200 0 120;         % 8. magenta
        40 150 200;        % 9. sky-blue
        0 160 160;         % 10 teal
        100 0 0;           % 11 dark red
        0 0 120;           % 12 navy
    ] / 255;

    if numCountries <= size(base,1)
        C = base(1:numCountries,:);
    else
        % if more colors needed → interpolate colormap
        extra = interpolate_colors(numCountries - size(base,1));
        C = [base; extra];
    end
end

function extra = interpolate_colors(k)
    % generate k extra distinct bright colors
    hsv = [linspace(0,1,k)' ones(k,1)*0.8 ones(k,1)*0.9];
    extra = hsv2rgb(hsv);
end
