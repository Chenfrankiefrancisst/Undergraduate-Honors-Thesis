function [HighD,LowD]=cumbootstrapcholcorrected(data,hor,c,iter,conf,p,n)

% Function to perform bootstrap method with Kilian (1998)'s bootstrap after
% bootstrap method:

% Author: Nicolo' Maffei Faccioli

% Step 1: 

[pi_hat,~,~,Y_initial,~,err]=VAR(data,p,c);

T=length(data)-length(Y_initial); % Data after losing lags

if c==1
    cons=pi_hat(1,:)';
    PI=pi_hat(2:end,:)';
else, cons=0;
    PI=pi_hat';
    
end

ExtraBigD=zeros(hor,n*n,2000);
Data_new=zeros(T+p,n,iter);
pi_hat_c=zeros(size(pi_hat,1),size(pi_hat,2));
err_c=zeros(T,n,iter);

for j=1:iter
    
% STEP 2: generate a new sample with bootstrap

Y_new=zeros(T,n);
Y_in= reshape(flipud(Y_initial)',1,[]);                  % 1*(n*p) row vector of [y_0 ... y_-p+1]

  for i=1:T
    Y_new(i,:)= cons' + Y_in*PI' + err(randi(T),:);  
    Y_in=[Y_new(i,:) Y_in(1:n*(p-1))];
  end
    Data_new(:,:,j)=[Y_initial ; Y_new];                        % Add the p initial lags to recreate the sample
  
% STEP 3 Correct for the biases:

[pi_hat_c(:,:,j),~,~,~,~,err_c(:,:,j)]=VAR(Data_new(:,:,j),p,c);

end

bias=mean(pi_hat_c,3)-pi_hat;
new_pi_hat=pi_hat-bias;

if c==1
    cons=new_pi_hat(1,:)';
    PI=new_pi_hat(2:end,:)';
else, cons=0;
    PI=new_pi_hat';  
end

[Y,X,Y_initial]=SUR(data,p,c);
Yfit=X*new_pi_hat; %Fitted value of Y
err=Y-Yfit; %Residuals

for j=1:2000
    
% STEP 2: generate a new sample with bootstrap

Y_new=zeros(T,n);
Y_in= reshape(flipud(Y_initial)',1,[]);                  % 1*(n*p) row vector of [y_0 ... y_-p+1]

  for i=1:T
    Y_new(i,:)= cons' + Y_in*PI' + err(randi(T),:);  
    Y_in=[Y_new(i,:) Y_in(1:n*(p-1))];
  end
    Data_new2=[Y_initial ; Y_new];                        % Add the p initial lags to recreate the sample
  
% STEP 3: Esimate VAR(p) with new sample and obtain IRF:

[pi_hat_c2,~,~,~,~,err_c2]=VAR(Data_new2,p,c);

if c==1
    BigA_c=[pi_hat_c2(2:end,:)'; eye(n*p-n) zeros(n*p-n,n)]; % BigA companion form, npxnp matrix
else, BigA_c=[pi_hat_c2'; eye(n*p-n) zeros(n*p-n,n)];
end

t=length(Data_new2)-n*p-n-p; % -p lags for n variables, -n constants, -p initally discarded observations
omega=(err_c2'*err_c2)./t; %estimate of omega
S=chol(omega,'lower'); %cholesky factorization, lower triangular matrix

C_c=zeros(n,n,hor);

for l=1:hor
    BigC_c=BigA_c^(l-1);
    C_c(:,:,l)=BigC_c(1:n,1:n)*S; % Impulse response functions of the Wold representation (Bootstrap)
end

C_c_1= reshape(permute(C_c,[3 2 1]),hor,n*n,[]);
ExtraBigD(:,:,j)=cumsum(C_c_1);

end

% Create bands:

LowD=zeros(hor,n*n);
HighD=zeros(hor,n*n);

for k=1:n
 for j=1:n
        Dmin = prctile(ExtraBigD(:,j+n*k-n,:),(100-conf)/2,3); %16th percentile
        LowD(:,j+n*k-n) = Dmin; %lower band
        Dmax = prctile(ExtraBigD(:,j+n*k-n,:),(100+conf)/2,3); %84th percentile
        HighD(:,j+n*k-n) = Dmax; %upper band
 end
end

end

