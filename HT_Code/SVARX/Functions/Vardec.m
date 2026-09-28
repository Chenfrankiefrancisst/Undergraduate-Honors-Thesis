err=finaldata(p+1:end,:)-X*pi_hat;
omega=cov(err); % Estimate of omega
S=chol(omega,'lower'); % Cholesky factorization, lower triangular matrix

D_var=zeros(n,n,hor);
D_ex=zeros(n,m,hor);

for j=1:hor
    
    BigD_var=(BigA^(j-1));
    D_var(:,:,j)=BigD_var(1:n,1:n)*S; % Impulse response functions of the Wold representation
    
    BigD_ex=BigD_var(1:n,1:n)*BigP(1:n,:);
    D_ex(:,:,j)=BigD_ex; % Impulse response functions of the Wold representation
end

D_wold_var=reshape(permute(D_var,[3 2 1]),hor,n*n,[]);
D_wold_ex=reshape(permute(D_ex,[3 2 1]),hor,n*m,[]);
D_mix=zeros(hor,n*n+n*m);

for i=1:n
D_mix(:,i*(n+m)-(n+m)+1:i*(n+m))=[D_wold_var(:,i*n-n+1:i*n) D_wold_ex(:,i*m-m+1:i*m)];
end

vardec=vdec_VARX(D_mix,hor,n,m,news);

% Plot:

vdec_plot_VARX;