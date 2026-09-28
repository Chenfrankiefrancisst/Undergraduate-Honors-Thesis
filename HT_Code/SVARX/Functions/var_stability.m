function [Acomp, SIGMA, U, V] = var_stability(Y, p)
% ---------------------------------------------------------
%   Computes VAR(p) stability using companion matrix
% ---------------------------------------------------------

[T, n] = size(Y);

% Build lag matrix
Ylag = [];
for i = 1:p
    Ylag = [Ylag, Y(p-i+1:end-i, :)];
end
Ynow = Y(p+1:end, :);

% Regression matrix
X = [ones(T-p,1), Ylag];

% OLS estimate
A = (X' * X) \ (X' * Ynow);    % (np+1)-by-n

% Strip constant
A1 = A(2:end, :)';            % n-by-(np)

% ---------------------------------------------------------
% Build companion matrix (np x np)
% ---------------------------------------------------------
Acomp = [A1 ;
         eye(n*(p-1)), zeros(n*(p-1), n)];

% Residuals
U = Ynow - X*A;
SIGMA = cov(U);

V = [];

end
