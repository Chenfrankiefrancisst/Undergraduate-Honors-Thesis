function [C,BigA,pi_hat,Y,X,Y_initial,Yfit,err] = VAR_irf(data,n,p,c,hor)

[pi_hat,Y,X,Y_initial,Yfit,err]=VAR(data,p,c); % VAR estimation

BigA=[pi_hat(2:end,:)'; eye(n*p-n) zeros(n*p-n,n)]; % np x np matrix

% Wold representation impulse responses:

C=zeros(n,n,hor);

for j=1:hor
    BigC=BigA^(j-1);
    C(:,:,j)=BigC(1:n,1:n); % IRF of the Wold representation
end

end