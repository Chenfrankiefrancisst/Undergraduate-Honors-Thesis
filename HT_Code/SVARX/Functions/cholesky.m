function [cirf]=cholesky(A,u,n,h,ind,stddata)

% Author: Sarah Zoi
% Inputs:
% - A: companion matrix
% - h: max horizon for impulse responses
% - n: number of endogenous variables
% - u: residual of the VAR model
% Output: 
% - nxn  matrix of cholesky irf
% - ind: index for the variables that need to cumulate 
% - stddata: vector of std of data as they are in the database(fd, sd,
% levels)
if isempty(stddata)==1 | nargin==5
    
    np=size(A,1);
    irf=nan(np,np,h);
    cirf=nan(n,n,h);
    CovU=cov(u);
    for i=1:h    
    irf(1:np,1:np,i)=A^(i-1);
    cirf(1:n,1:n,i)=irf(1:n,1:n,i)*chol(CovU,'lower');
    end
    % wold IRFs
    irf=irf(1:n,1:n,:);
    % cumulate cholesky IRFs
    cirf(ind,:,:)=cumsum(cirf(ind,:,:),3);

    else
            stddata=repmat(stddata',1,n);
    np=size(A,1);
    irf=nan(np,np,h);
    cirf=nan(n,n,h);
    CovU=cov(u);
    for i=1:h    
    irf(1:np,1:np,i)=A^(i-1);
    cirf(1:n,1:n,i)=irf(1:n,1:n,i)*chol(CovU,'lower').*stddata;
    end
    % wold IRFs
    irf=irf(1:n,1:n,:);
    % cumulate cholesky IRFs
    cirf(ind,:,:)=cumsum(cirf(ind,:,:),3);      
end
end