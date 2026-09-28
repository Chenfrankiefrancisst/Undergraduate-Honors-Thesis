function [beta, se, tstat, pval] = newey_west(X, y, L)

    % OLS
    beta = (X' * X) \ (X' * y);
    T = size(X,1);
    u = y - X * beta;
    
    % Meat
    S = zeros(size(X,2));
    for t = 1:T
        S = S + (X(t,:)' * X(t,:)) * (u(t)^2);
    end

    % Add lag terms
    for l = 1:L
        w = 1 - l/(L+1);
        Gamma = zeros(size(X,2));
        for t = l+1:T
            Gamma = Gamma + (X(t,:)' * X(t-l,:)) * (u(t) * u(t-l));
        end
        S = S + w * (Gamma + Gamma');
    end
    
    % Covariance matrix
    V = (X' * X) \ S / (X' * X);

    % Standard errors
    se = sqrt(diag(V));
    
    % Test statistics
    tstat = beta ./ se;
    pval = 2*(1 - tcdf(abs(tstat), T - size(X,2)));
end
