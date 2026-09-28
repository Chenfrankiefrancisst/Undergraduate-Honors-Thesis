if q<=p
    newsX=[news(q+1:end,:),lagmaker(news,q)];
    newsX=newsX(p-q+1:end,:);
    Y=finaldata(p+1:end,:);
else
    newsX=[news(q+1:end,:),lagmaker(news,q)];
    Y=finaldata(q+1:end,:);
end

X=[ones(length(lagmaker(finaldata,p)),1),lagmaker(finaldata,p),newsX];

if restr==1

X_1=X(:,1:end-m);
X_2=X;

pi_hat_1=(X_1'*X_1)\X_1'*Y(:,1:pos_shock-1);
pi_hat_1=[pi_hat_1;zeros(m,pos_shock-1)];
pi_hat_2=(X_2'*X_2)\X_2'*Y(:,pos_shock:end);
pi_hat=[pi_hat_1,pi_hat_2];

elseif restr==2
    
X_1=X(:,[1:end-m end-m+2:end]);
X_2=X;

pi_hat_1=(X_1'*X_1)\X_1'*Y(:,1:pos_shock-1);
pi_hat_1=[pi_hat_1(1:end-m+1,:);zeros(1,pos_shock-1);pi_hat_1(end-m+2:end,:)];
pi_hat_2=(X_2'*X_2)\X_2'*Y(:,pos_shock:end);
pi_hat=[pi_hat_1,pi_hat_2];

elseif restr==3 

X_1=X(:,1:end-m);
X_2=X(:,1:end-m+1);
X_3=X;

pi_hat_1=(X_1'*X_1)\X_1'*Y(:,1:pos_shock-1);
pi_hat_1=[pi_hat_1;zeros(m,pos_shock-1)];
pi_hat_2=(X_2'*X_2)\X_2'*Y(:,pos_shock);
pi_hat_2=[pi_hat_2;zeros(m-1,1)];
pi_hat_3=(X_3'*X_3)\X_3'*Y(:,pos_shock+1:end);
pi_hat=[pi_hat_1,pi_hat_2,pi_hat_3];

elseif restr==4
    
X_1=X(:,[1:end-m end-m+2:end]);
X_2=X(:,1:end-m+1);
X_3=X;

pi_hat_1=(X_1'*X_1)\X_1'*Y(:,1:pos_shock-1);
pi_hat_1=[pi_hat_1(1:end-m+1,:);zeros(1,pos_shock-1);pi_hat_1(end-m+2:end,:)];
pi_hat_2=(X_2'*X_2)\X_2'*Y(:,pos_shock);
pi_hat_2=[pi_hat_2;zeros(1,1)];
pi_hat_3=(X_3'*X_3)\X_3'*Y(:,pos_shock+1:end);
pi_hat=[pi_hat_1,pi_hat_2,pi_hat_3];

else
    
pi_hat=(X'*X)\X'*Y;    
    
end

BigA=[pi_hat(2:end-m*(q+1),:)'; eye(n*p-n) zeros(n*p-n,n)]; % BigA companion form, npxnp matrix
BigP=[pi_hat(end-m*(q+1)+1:end,:)';zeros(n*p-n,m*(q+1))];

% IRFs:

C=zeros(n,m,hor);

for j=1:hor
    BigC=(BigA^(j-1))*BigP;
    C(:,:,j)=BigC(1:n,1:m); 
end

C_wold=reshape(permute(C,[3 2 1]),hor,n*m,[]);

for i=[1,2,6]
C_wold(:,i*m-1:i*m)=cumsum(C_wold(:,i*m-1:i*m),1);
end

% Bootstrap:

iter=1000;
conf=68;
trend=0;

[HighC,LowC,ExtraBigC]=bootstrapSVARX_diff(finaldata,news,hor,1,trend,iter,conf,p,q,n,cum,pos_shock,restr,nonlin); % function that performs bootstrap
[HighC90,LowC90]=bootstrapSVARX_diff(finaldata,news,hor,1,trend,iter,90,p,q,n,cum,pos_shock,restr,nonlin); % function that performs bootstrap