function [HighC,LowC,HighC_magnitude,LowC_magnitude]=bootstrapSVARX_pos(data,news,hor,c,trend,iter,conf,p,q,n,magnitude,cum,pos_shock,restr,nonlin)

% Function to perform bootstrap method
% Author: Nicolo' Maffei Faccioli

% Step 1: 

Y_initial=data(1:p,:);
m=size(news,2);

if q<=p
    newsX_nico=[news(q+1:end,:),lagmaker(news,q)];
    newsX_nico=newsX_nico(p-q+1:end,:);
    Y_nico=data(p+1:end,:);
else 
    newsX_nico=[news(q+1:end,:),lagmaker(news,q)];
    Y_nico=data(q+1:end,:);
end

if trend==1
X=[ones(size(lagmaker(data,p),1),1),cumsum(ones(size(lagmaker(data,p),1),1)).^2,lagmaker(data,p),newsX_nico];
else, X=[ones(size(lagmaker(data,p),1),1),lagmaker(data,p),newsX_nico];
end

if restr==1

X_1=X(:,1:end-m);
X_2=X;

pi_hat_1=(X_1'*X_1)\X_1'*Y_nico(:,1:pos_shock-1);
pi_hat_1=[pi_hat_1;zeros(m,pos_shock-1)];
pi_hat_2=(X_2'*X_2)\X_2'*Y_nico(:,pos_shock:end);
pi_hat=[pi_hat_1,pi_hat_2];

elseif restr==2 % only shock restricted

X_1=X(:,[1:end-m end-m+2:end]);
X_2=X;

pi_hat_1=(X_1'*X_1)\X_1'*Y_nico(:,1:pos_shock-1);
pi_hat_1=[pi_hat_1(1:end-m+1,:);zeros(1,pos_shock-1);pi_hat_1(end-m+2:end,:)];
pi_hat_2=(X_2'*X_2)\X_2'*Y_nico(:,pos_shock:end);
pi_hat=[pi_hat_1,pi_hat_2];

elseif restr==3 % both shock and square restricted for EBP only

X_1=X(:,1:end-m);
X_2=X(:,1:end-m+1);
X_3=X;

pi_hat_1=(X_1'*X_1)\X_1'*Y_nico(:,1:pos_shock-1);
pi_hat_1=[pi_hat_1;zeros(m,pos_shock-1)];
pi_hat_2=(X_2'*X_2)\X_2'*Y_nico(:,pos_shock);
pi_hat_2=[pi_hat_2;zeros(1,1)];
pi_hat_3=(X_3'*X_3)\X_3'*Y_nico(:,pos_shock+1:end);
pi_hat=[pi_hat_1,pi_hat_2,pi_hat_3];

elseif restr==4
    
X_1=X(:,[1:end-m end-m+2:end]);
X_2=X(:,1:end-m+1);
X_3=X;

pi_hat_1=(X_1'*X_1)\X_1'*Y_nico(:,1:pos_shock-1);
pi_hat_1=[pi_hat_1(1:end-m+1,:);zeros(1,pos_shock-1);pi_hat_1(end-m+2:end,:)];
pi_hat_2=(X_2'*X_2)\X_2'*Y_nico(:,pos_shock);
pi_hat_2=[pi_hat_2;zeros(1,1)];
pi_hat_3=(X_3'*X_3)\X_3'*Y_nico(:,pos_shock+1:end);
pi_hat=[pi_hat_1,pi_hat_2,pi_hat_3];

else
    
pi_hat=(X'*X)\X'*data(p+1:end,:);

end

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
Y_in=reshape(flipud(Y_initial)',1,[]) ; % 1*(n*p) row vector of [y_0 ... y_-p+1]
newsX_new=zeros(T,m);

  for i=1:T
    rnd=randi(T);
    Y_new(i,:)=((length(Y_initial)+i).^2).*(t)' + cons' + Y_in*PI' + newsX_nico(rnd,:)*newsPI' + err(rnd,:);  
    Y_in=[Y_new(i,:) Y_in(1:n*(p-1))];
    newsX_new(i,:)=newsX_nico(rnd,1:m);
  end
  
Data_new=[Y_initial ; Y_new];   % Add the p initial lags to recreate the sample
    
% STEP 3: Esimate VAR and VARX with new sample and obtain IRF:

[~,S_new,~,~,~,~,~,~,~,err_new]=chol_irf(Data_new,n,p,c,hor); % estimate the first step

eta_new=(S_new\err_new')'; % Cholesky shocks

if nonlin==1
news_new=[news(1:p,:); eta_new(:,pos_shock) eta_new(:,pos_shock).^2];
else, news_new=[news(1:p,:); eta_new(:,pos_shock) abs(eta_new(:,pos_shock))];
end

if q<=p
newsX=[news_new(q+1:end,:),lagmaker(news_new,q)];
newsX=newsX(p-q+1:end,:);
else, newsX=[news_new(q+1:end,:),lagmaker(news_new,q)];
end

if trend==1
X_c=[ones(size(lagmaker(Data_new,p),1),1),cumsum(ones(size(lagmaker(Data_new,p),1),1)).^2,lagmaker(Data_new,p),newsX];
else, X_c=[ones(size(lagmaker(Data_new,p),1),1),lagmaker(Data_new,p),newsX];
end

if restr==1

X_1_c=X_c(:,1:end-m);
X_2_c=X_c;

pi_hat_1_c=(X_1_c'*X_1_c)\X_1_c'*Data_new(p+1:end,1:pos_shock-1);
pi_hat_1_c=[pi_hat_1_c;zeros(m,pos_shock-1)];
pi_hat_2_c=(X_2_c'*X_2_c)\X_2_c'*Data_new(p+1:end,pos_shock:end);
pi_hat_c=[pi_hat_1_c,pi_hat_2_c];

elseif restr==2

X_1_c=X_c(:,[1:end-m end-m+2:end]);
X_2_c=X_c;

pi_hat_1_c=(X_1_c'*X_1_c)\X_1_c'*Data_new(p+1:end,1:pos_shock-1);
pi_hat_1_c=[pi_hat_1_c(1:end-m+1,:);zeros(1,pos_shock-1);pi_hat_1_c(end-m+2:end,:)];
pi_hat_2_c=(X_2_c'*X_2_c)\X_2_c'*Data_new(p+1:end,pos_shock:end);
pi_hat_c=[pi_hat_1_c,pi_hat_2_c];

elseif restr==3

X_1_c=X_c(:,1:end-m);
X_2_c=X_c(:,1:end-m+1);
X_3_c=X_c;

pi_hat_1_c=(X_1_c'*X_1_c)\X_1_c'*Data_new(p+1:end,1:pos_shock-1);
pi_hat_1_c=[pi_hat_1_c;zeros(m,pos_shock-1)];
pi_hat_2_c=(X_2_c'*X_2_c)\X_2_c'*Data_new(p+1:end,pos_shock);
pi_hat_2_c=[pi_hat_2_c;zeros(1,1)];
pi_hat_3_c=(X_3_c'*X_3_c)\X_3_c'*Data_new(p+1:end,pos_shock+1:end);
pi_hat_c=[pi_hat_1_c,pi_hat_2_c,pi_hat_3_c];

elseif restr==4

X_1_c=X_c(:,[1:end-m end-m+2:end]);
X_2_c=X_c(:,1:end-m+1);
X_3_c=X_c;

pi_hat_1_c=(X_1_c'*X_1_c)\X_1_c'*Data_new(p+1:end,1:pos_shock-1);
pi_hat_1_c=[pi_hat_1_c(1:end-m+1,:);zeros(1,pos_shock-1);pi_hat_1_c(end-m+2:end,:)];
pi_hat_2_c=(X_2_c'*X_2_c)\X_2_c'*Data_new(p+1:end,pos_shock);
pi_hat_2_c=[pi_hat_2_c;zeros(1,1)];
pi_hat_3_c=(X_3_c'*X_3_c)\X_3_c'*Data_new(p+1:end,pos_shock+1:end);
pi_hat_c=[pi_hat_1_c,pi_hat_2_c,pi_hat_3_c];

else

pi_hat_c=(X_c'*X_c)\X_c'*Data_new(p+1:end,:);

end

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

ExtraBigC_magnitude(:,:,j)=magnitude.*C_c_1(:,1:m:end)+(magnitude.^2).*C_c_1(:,2:m:end);

end

% Create bands:

LowC=zeros(hor,n*m);
HighC=zeros(hor,n*m);

LowC_magnitude=zeros(hor,n);
HighC_magnitude=zeros(hor,n);

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
 
if cum==0

 for j=1:n
        Cmin = prctile(ExtraBigC_magnitude(:,j,:),(100-conf)/2,3); %16th percentile
        LowC_magnitude(:,j) = Cmin; %lower band
        Cmax = prctile(ExtraBigC_magnitude(:,j,:),(100+conf)/2,3); %84th percentile
        HighC_magnitude(:,j) = Cmax; %upper band
 end
 
elseif cum==1
    
 for j=1:n
        Cmin = prctile(cumsum(ExtraBigC_magnitude(:,j,:),1),(100-conf)/2,3); %16th percentile
        LowC_magnitude(:,j) = Cmin; %lower band
        Cmax = prctile(cumsum(ExtraBigC_magnitude(:,j,:),1),(100+conf)/2,3); %84th percentile
        HighC_magnitude(:,j) = Cmax; %upper band
 end
    
end

end

