function [HighD,LowD]=chol_boot(data,hor,c,iter,conf,p,n,cum)

% Function to perform bootstrap method with Kilian (1998)'s bootstrap after
% bootstrap method:

% Author: Nicolo' Maffei Faccioli

% Step 1: Estimate the VAR(p) and store the residuals

[pi_hat,~,~,Y_initial,~,err]=VAR(data,p,c);

if c==1
    cons=pi_hat(1,:)';
    PI=pi_hat(2:end,:)';
else, cons=0;
    PI=pi_hat';    
end

% Set up matrices for bootstrap loop:

T=length(data)-p; % Data after losing lags

ExtraBigD=zeros(hor,n*n,2000);
Data_new=zeros(T+p,n,iter);
pi_hat_c=zeros(size(pi_hat,1),size(pi_hat,2));
err_c=zeros(T,n,iter);

for j=1:iter
    
% STEP 2: Generate a new sample for each j with bootstrap:

Y_new=zeros(T,n);
Y_in=reshape(flipud(Y_initial)',1,[]);  % 1 x n*p row vector of [y_0,...,y_-p+1]

  for i=1:T % create the new sample 
    Y_new(i,:)= cons' + Y_in*PI' + err(randi(T),:);  
    Y_in=[Y_new(i,:) Y_in(1:n*(p-1))];
  end
  
    Data_new(:,:,j)=[Y_initial ; Y_new]; % Add the p initial lags to recreate the sample
  
% STEP 3: Estimate VAR(p) with the newly created j=1,...,iter samples:

[pi_hat_c(:,:,j),~,~,~,~,err_c(:,:,j)]=VAR(Data_new(:,:,j),p,c);

end

% STEP 4: Correct for the bias in the IRFs (Kilian, 1998):

bias=mean(pi_hat_c,3)-pi_hat;
new_pi_hat=pi_hat-bias;

if c==1
    cons=new_pi_hat(1,:)';
    PI=new_pi_hat(2:end,:)';
else, cons=0;
    PI=new_pi_hat';  
end

% STEP 5: Perform the bootstrap procedure again, with the bias correction:

[Y,X,Y_initial]=SUR(data,p,c);
Yfit=X*new_pi_hat; %Fitted value of Y
err=Y-Yfit; %Residuals

for j=1:2000
    
% Generate a new sample with bootstrap

Y_new=zeros(T,n);
Y_in=reshape(flipud(Y_initial)',1,[]); % 1 x n*p row vector of [y_0,...,y_-p+1]

  for i=1:T
    Y_new(i,:)= cons' + Y_in*PI' + err(randi(T),:);  
    Y_in=[Y_new(i,:) Y_in(1:n*(p-1))];
  end
    Data_new2=[Y_initial ; Y_new];  % Add the p initial lags to recreate the sample
  
% Esimate VAR(p) with new samples and obtain Cholesky IRFs:

D=chol_irf(Data_new2,n,p,c,hor);

MiddleD=(reshape(permute(D,[3 2 1]),hor,n*n,[]));

ExtraBigD(:,:,j)=MiddleD;

end

% Create bands:

LowD=zeros(hor,n*n);
HighD=zeros(hor,n*n);

if cum==1

for k=1:n
 for j=1:n
        Dmin = prctile(cumsum(ExtraBigD(:,j+n*k-n,:),1),(100-conf)/2,3); %16th percentile
        LowD(:,j+n*k-n) = Dmin; %lower band
        Dmax = prctile(cumsum(ExtraBigD(:,j+n*k-n,:),1),(100+conf)/2,3); %84th percentile
        HighD(:,j+n*k-n) = Dmax; %upper band
 end
end

else
    
for k=1:n
 for j=1:n
        Dmin = prctile(ExtraBigD(:,j+n*k-n,:),(100-conf)/2,3); %16th percentile
        LowD(:,j+n*k-n) = Dmin; %lower band
        Dmax = prctile(ExtraBigD(:,j+n*k-n,:),(100+conf)/2,3); %84th percentile
        HighD(:,j+n*k-n) = Dmax; %upper band
 end
end  
    
end

end

