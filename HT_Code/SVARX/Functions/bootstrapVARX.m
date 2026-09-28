function [HighC,LowC]=bootstrapVARX(data,news,hor,c,trend,iter,conf,p,q,n,cum)

% Function to perform bootstrap method
% Author: Nicolo' Maffei Faccioli

% Step 1: 

Y_initial=data(1:p,:);
m=size(news,2);

if q<=p
newsX=[news(q+1:end,:),lagmaker(news,q)];
newsX=newsX(p-q+1:end,:);
else, newsX=[news(q+1:end,:),lagmaker(news,q)];
end

if trend==1
X=[ones(size(lagmaker(data,p),1),1),cumsum(ones(size(lagmaker(data,p),1),1)).^2,lagmaker(data,p),newsX];

else, X=[ones(size(lagmaker(data,p),1),1),lagmaker(data,p),newsX];

end

pi_hat=(X'*X)\X'*data(p+1:end,:);
err=data(p+1:end,:)-X*pi_hat;

T=length(data)-size(Y_initial,1); % Data after losing lags

if c==1
    cons=pi_hat(1,:)';
    PI=pi_hat(2+trend:end,:)';
else, cons=0;
    PI=pi_hat(1+trend,:)';   
end

if trend==1
    t=pi_hat(1+c,:)';
else, t=0;
end

PI=PI';
newsPI=PI(end-m*(q+1)+1:end,:)';
PI=PI(1:end-m*(q+1),:)';

ExtraBigC=zeros(hor,n*m,iter);

for j=1:iter
    
% STEP 2: generate a new sample with bootstrap

Y_new=zeros(T,n);
Y_in= reshape(flipud(Y_initial)',1,[]) ; % 1*(n*p) row vector of [y_0 ... y_-p+1]

  for i=1:T
    Y_new(i,:)= ((length(Y_initial)+i).^2).*(t)' + cons' + Y_in*PI' + newsX(i,:)*newsPI' + err(randi(T),:);  
    Y_in=[Y_new(i,:) Y_in(1:n*(p-1))];
  end
  
    Data_new=[Y_initial ; Y_new];                        % Add the p initial lags to recreate the sample
    
% STEP 3: Esimate VAR(p) with new sample and obtain IRF:

if q<=p
newsX=[news(q+1:end,:),lagmaker(news,q)];
newsX=newsX(p-q+1:end,:);
else, newsX=[news(q+1:end,:),lagmaker(news,q)];
end

if trend==1
X_c=[ones(size(lagmaker(Data_new,p),1),1),cumsum(ones(size(lagmaker(Data_new,p),1),1)).^2,lagmaker(Data_new,p),newsX];

else, X_c=[ones(size(lagmaker(Data_new,p),1),1),lagmaker(Data_new,p),newsX];

end

pi_hat_c=(X_c'*X_c)\X_c'*Data_new(p+1:end,:);

if c==1
BigA_c=[pi_hat_c(2+trend:end-m*(q+1),:)'; eye(n*p-n) zeros(n*p-n,n)]; % BigA companion form, np+1 x np+1 matrix
else, BigA_c=[pi_hat_c(1+trend:end-m*(q+1),:)'; eye(n*p-n) zeros(n*p-n,n)]; % BigA companion form, np x np matrix
end
BigP=[pi_hat_c(end-m*(q+1)+1:end,:)';zeros(n*p-n,m*(q+1))];

C_c=zeros(n,m,hor);

for l=1:hor
    BigC_c=BigA_c^(l-1)*BigP;
    C_c(:,:,l)=BigC_c(1:n,1:m); % Impulse response functions of the Wold representation (Bootstrap)
end

C_c_1= reshape(permute(C_c,[3 2 1]),hor,n*m,[]);
ExtraBigC(:,:,j)=C_c_1;

end

% Create bands:

LowC=zeros(hor,n*m);
HighC=zeros(hor,n*m);

if cum==0

 for j=1:n*m
        Cmin = prctile(ExtraBigC(:,j,:),(100-conf)/2,3); %16th percentile
        LowC(:,j) = Cmin; %lower band
        Cmax = prctile(ExtraBigC(:,j,:),(100+conf)/2,3); %84th percentile
        HighC(:,j) = Cmax; %upper band
 end
 
elseif cum==1
    
 for j=1:n*m
        Cmin = prctile(cumsum(ExtraBigC(:,j,:),1),(100-conf)/2,3); %16th percentile
        LowC(:,j) = Cmin; %lower band
        Cmax = prctile(cumsum(ExtraBigC(:,j,:),1),(100+conf)/2,3); %84th percentile
        HighC(:,j) = Cmax; %upper band
 end
    
end

end

