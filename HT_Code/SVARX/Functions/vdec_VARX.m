function vardec = vdec_VARX(D,hor,n,m,news)

MiddleDsquare=D.^2;

for i=n+m:n+m:n*n+n*m
MiddleDsquare(:,i)=MiddleDsquare(:,i).*var(news(:,m));
end


denom=zeros(hor,n);

for k=[1:n;1:n+m:n*n+n*m]
    
    denom(:,k(1))=cumsum(sum(MiddleDsquare(:,k(2):k(2)+n+m-1),2));
    
end

denomtot=zeros(hor,n*n+n*m);

for k=[1:n;1:n+m:n*n+n*m]
    
    denomtot(:,k(2):k(2)+n+m-1)=denom(:,k(1)).*ones(hor,n+m);
    
end

vardec=zeros(hor,n*n+n*m);

for j=1:n*n+n*m
    
    vardec(:,j)=cumsum(MiddleDsquare(:,j))./denomtot(:,j);
    
end
