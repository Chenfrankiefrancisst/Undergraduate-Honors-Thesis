% Compute the historical decomposition:

err=finaldata(p+1:end,:)-X*pi_hat; % reduced-form errors of VARX
omega=cov(err); % Estimate of omega
S=chol(omega,'lower'); % Cholesky factorization, lower triangular matrix

D_var=zeros(n,n,hor); % IRFs othogonal shocks endogenous variables
D_ex=zeros(n,m,hor); % IRFs exogenous variables

for j=1:hor
    
    BigD_var=(BigA^(j-1));
    D_var(:,:,j)=BigD_var(1:n,1:n)*S;
    
    BigD_ex=BigD_var(1:n,1:n)*BigP(1:n,:);
    D_ex(:,:,j)=BigD_ex; 
end

D_wold_var=reshape(permute(D_var,[3 2 1]),hor,n*n,[]);
D_wold_ex=reshape(permute(D_ex,[3 2 1]),hor,n*m,[]);

% Historical decomposition:

gamma=(S\err')'; % orthogonal shocks VARX

hist_var=zeros(size(gamma,1),n,n);

for i=1:n
    
    for j=1:n
    
    hist_var(:,i,j)=filter(squeeze(D_var(i,j,:)),1,gamma(:,j));
    
    end
    
end


hist_l=zeros(size(news(p+1:end,:),1),n); % Linear term financial shock
hist_nl=zeros(size(news(p+1:end,:),1),n); % Nonlinear term financial shock

for i=1:n

hist_l(:,i)=filter(squeeze(D_ex(i,1,:)),1,news(p+1:end,1));
hist_nl(:,i)=filter(squeeze(D_ex(i,2,:)),1,news(p+1:end,2));

end

shock_contr=squeeze(sum(hist_var,3))+hist_l+hist_nl; % stochastic component of VARX, i.e. part explained by shocks

totdec=[squeeze(sum(hist_var,3)) hist_l hist_nl];

positive=totdec>=0;
negative=totdec<0;
histdecpos=zeros(size(totdec,1),size(totdec,2));
histdecneg=zeros(size(totdec,1),size(totdec,2));

for t=1:size(totdec,1)
    
    for k=1:size(totdec,2)
    
    if positive(t,k)==1
        
        histdecpos(t,k)=totdec(t,k); 
        
    else, histdecpos(t,k)=NaN;
        
    end
    
    end
    
end

for t=1:size(totdec,1)
    
    for k=1:size(totdec,2)
    
    if positive(t,k)==0
        
        histdecneg(t,k)=totdec(t,k); 
        
    else, histdecneg(t,k)=NaN;
        
    end
    
    end
    
end