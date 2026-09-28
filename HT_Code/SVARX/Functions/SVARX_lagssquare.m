if squareonly==1
    
if q<=p
    newsX=[news(q+1:end,:),lagmaker(news(:,2),q)];
    newsX=newsX(p-q+1:end,:);
    Y=finaldata(p+1:end,:);
else
    newsX=[news(q+1:end,:),lagmaker(news(:,2),q)];
    Y=finaldata(q+1:end,:);
end
    
    
else

if q<=p
    newsX=[news(q+1:end,:),lagmaker(news,q)];
    newsX=newsX(p-q+1:end,:);
    Y=finaldata(p+1:end,:);
else
    newsX=[news(q+1:end,:),lagmaker(news,q)];
    Y=finaldata(q+1:end,:);
end

end

X=[ones(length(lagmaker(finaldata,p)),1),lagmaker(finaldata,p),newsX];

if restr==1

X_1=X(:,[1:n*p+1 n*p+2+m:end]);
X_2=X;

pi_hat_1=(X_1'*X_1)\X_1'*Y(:,1:pos_shock-1);
%pi_hat_1=[pi_hat_1;zeros(m,pos_shock-1)];
pi_hat_1=[pi_hat_1(1:n*p+1,:);zeros(m,pos_shock-1);pi_hat_1(n*p+2:end,:)];
pi_hat_2=(X_2'*X_2)\X_2'*Y(:,pos_shock:end);
pi_hat=[pi_hat_1,pi_hat_2];

elseif restr==2

X_1=X(:,[1:n*p+1 n*p+1+m:end]);
X_2=X;

pi_hat_1=(X_1'*X_1)\X_1'*Y(:,1:pos_shock-1);
%pi_hat_1=[pi_hat_1;zeros(m,pos_shock-1)];
pi_hat_1=[pi_hat_1(1:n*p+1,:);zeros(1,pos_shock-1);pi_hat_1(n*p+2:end,:)];
pi_hat_2=(X_2'*X_2)\X_2'*Y(:,pos_shock:end);
pi_hat=[pi_hat_1,pi_hat_2];

elseif restr==3

X_1=X(:,[1:n*p+1 n*p+2+m:end]);
X_2=X(:,[1:n*p+2 n*p+2+m:end]);
X_3=X;

pi_hat_1=(X_1'*X_1)\X_1'*Y(:,1:pos_shock-1);
%pi_hat_1=[pi_hat_1;zeros(m,pos_shock-1)];
pi_hat_1=[pi_hat_1(1:n*p+1,:);zeros(m,pos_shock-1);pi_hat_1(n*p+2:end,:)];
pi_hat_2=(X_2'*X_2)\X_2'*Y(:,pos_shock);
pi_hat_2=[pi_hat_2(1:n*p+2,:);zeros(m-1,1);pi_hat_2(n*p+2+1:end,:)];
pi_hat_3=(X_3'*X_3)\X_3'*Y(:,pos_shock+1:end);
pi_hat=[pi_hat_1,pi_hat_2,pi_hat_3];

else

pi_hat=(X'*X)\X'*Y;       
    
end

if squareonly==1
    
BigA=[pi_hat(2:end-(m-1)*(q+1)-1,:)'; eye(n*p-n) zeros(n*p-n,n)]; % BigA companion form, npxnp matrix
BigP=[pi_hat(end-(m-1)*(q+1):end,:)';zeros(n*p-n,(m-1)*(q+1)+1)];

else

BigA=[pi_hat(2:end-m*(q+1),:)'; eye(n*p-n) zeros(n*p-n,n)]; % BigA companion form, npxnp matrix
BigP=[pi_hat(end-m*(q+1)+1:end,:)';zeros(n*p-n,m*(q+1))];

end

% Wold representation impulse responses:

C=zeros(n,m,hor);

for j=1:hor
    BigC=(BigA^(j-1))*BigP;
    C(:,:,j)=BigC(1:n,1:m); % Impulse response functions of the Wold representation
end

C_wold=reshape(permute(C,[3 2 1]),hor,n*m,[]);

if cum==1
C_wold=cumsum(C_wold,1);
end

% Bootstrap:

iter=1000;
conf=68;

trend=0;

[HighC,LowC]=bootstrapVARX_new_lagssquare(finaldata,news,hor,1,trend,iter,conf,p,q,n,cum,pos_shock,restr,squareonly); % function that performs bootstrap
[HighC90,LowC90]=bootstrapVARX_new_lagssquare(finaldata,news,hor,1,trend,iter,90,p,q,n,cum,pos_shock,restr,squareonly); % function that performs bootstrap