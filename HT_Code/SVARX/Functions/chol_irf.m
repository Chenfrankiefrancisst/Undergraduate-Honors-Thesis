function [D,S,C,BigA,pi_hat,Y,X,Y_initial,Yfit,err,eta] = chol_irf(data,n,p,c,hor)

[pi_hat,Y,X,Y_initial,Yfit,err]=VAR(data,p,c); % VAR estimation

% Companion matrix representation:

BigA=[pi_hat(2:end,:)'; eye(n*p-n) zeros(n*p-n,n)]; % np x np matrix

% Wold representation IRFs:

C=zeros(n,n,hor); % reduced-form IRFs

for j=1:hor
    BigC=BigA^(j-1);
    C(:,:,j)=BigC(1:n,1:n); % IRFs of the Wold representation
end

% Cholesky IRFs:

omega=cov(err); % Estimate of omega
S=chol(omega,'lower'); % Cholesky factorization, lower triangular matrix

D=zeros(n,n,hor); % Cholesky IRFs

for j=1:hor
    D(:,:,j)=C(:,:,j)*S; 
end

eta=(S\err')'; % Cholesky shocks

end