function [cirfboot,cirf, Ac, lb, ub]=cholboot(data,p,const,h,ind,maxboot,cb,gr,stddata)
% Author: Sarah Zoi

% Note: this function works together with functions'VAR.m' and 'cholesky.m'

% Inputs:
% - data: columns are different variables (T*n)
% - p: n of lags
% - const: 1 if a constant is required in the model
% - h: horizon for impulse responses
% - ind: indices for those variables in levels for which IRFs need to be
%       cumulated

% - maxboot: maximum number of repetitions for the bootstrap
% - gr: 1 if graph is required

% Output:
% - cirfboot: cholesky impulse response functions maxboot cubic matrices
% of dimension (n*n*h)
% - cirf: cholesky impulse response
% - Ac: matrix of VAR coefficients
% - lb, ub: lower and upper bounds for confidence interval
if nargin==8 & isempty(ind)==0
   stddata=[];
elseif nargin==8 & isempty(ind)==1
    ind=[];
end
n= size(data,2);
T=size(data,1);
% estimate the VAR the first time and compute residuals
[X,Y, PI, Ac, u] = VAR(data,p, const);

% compute IRFs
[cirf]=cholesky(Ac,u,n,h,ind,stddata);
k=maxboot;
% Bootstrapping
cirfboot=nan(n,n,h,k);

for i=1:k
    
    % generate space for new data series
    Yfake=nan(T-p,n); 
    
    if const==1
   Xfake = ones(T-p,n*p+1);
    else
        Xfake = ones(T-p,n*p);
    end

        % generate new data series
        for j=1:T-p
        index=randi([1, T-p]);
          if j==1
          Xfake(j,:)=X(1,:);
          else
            Xfake(j,:)=[ones(const,1) Yfake(j-1,:) Xfake(j-1,1+const:n*(p-1)+const)];
          end
        Yfake(j,:)=Xfake(j,:)*PI+u(index,:);
        end
        datafake=Yfake;
    
    % run the VAR again and find cholesky impulse responses
    [Xtemp,Ytemp, PItemp, Atemp, utemp] = VAR(datafake,p, const);
    [cirftemp]=cholesky(Atemp,utemp,n,h,ind,stddata);
    % store the irf for iteration i
    cirfboot(1:n,1:n,:,i)=cirftemp;
end
% compute percentiles
% generate space to store percentiles
lb=nan(n,n,h);
ub=nan(n,n,h);
meanirf=nan(n,n,h);
for i=1:n
    for j=1:n
lb(i,j,:)=prctile(cirfboot(i,j,:,:),cb(1),4);
ub(i,j,:)=prctile(cirfboot(i,j,:,:),cb(2),4);
meanirf(i,j,:)=prctile(cirfboot(i,j,:,:),50,4);
    end
end


% plot if required
if gr==1
    figure
    g=0;
    for i=1:n  
        for j=1:n
            g=g+1;
    subplot(n,n,g), 
    plot(squeeze(cirf(i,j,:)),'k')
    %plot(squeeze(prctile(cirfboot(i,j,:,:),50,4)),'k')
    hold on
    plot(squeeze(lb(i,j,:)),':k')
    hold on
    plot(squeeze(ub(i,j,:)),':k')
    %plot(squeeze(prctile(cirfboot(i,j,:,:),cb(2),4)),':k')
    grid on
    axis tight
    %title([names{j},' on ',names{i}])
        end
    end
    
else
    
end

end