function [aicL, bicL, hqcL] = AicBicHqc(Data, pmax, varargin)
% The function [aicL, bicL, hqcL] = AicBicHqc(Data, pmax, varargin) gives 
% the AIC and BIC and HQC information criteria for lag selection.
%
% AIC=logV+2d/T 
% BIC=logV+logT*d/T
% HQC=logV+2*d/T*log(log(T));
%
% where:
% V=det(1/T Sum_{t=1}^T e_t*e_t')
% d # parameters 
% T # of observations
%
% INPUTS:
% Data - T*N data sample
% pmax - scalar number of lags to be checked
% varargin{1} - string 'tableON' to get the table (default OFF)
%
% OUTPUTS:
% aicL - number of lags suggested by aic
% bicL - number of lags suggested by bic
% hqcL - number of lags suggested by bic

% Check input:
narginchk(2,3)

if nargin==3 && ischar(varargin{1});
    Table=varargin{1};
else
    Table='tableOFF';
end

% Read Data:
[T,N] = size(Data);
constant=1;
% Preallocate for looping:
aic=NaN(pmax,1);
bic=NaN(pmax,1);
hqc=NaN(pmax,1);
for p = 1:pmax
    [~,~,~,~,~,e]=varestimate(Data,p,constant);
    V = det((e'*e)/T);  % Estimated varcov
    %d = p*N^2;          % Number of parameters
    d = p*N^2 + N;      % Number of parameters
    aic(p) = log(V) + 2*d/T;
    bic(p) = log(V) + log(T)*d/T;
    hqc(p) = log(V) + 2*d/T*log(log(T)); 
    % Check that the varcov is positive definite
    [~,test]=chol(e'*e);
    if test~=0 % if not discard the lag option
        aic(p)=Inf; bic(p)=Inf; hqc(p)=Inf;
    end
end
[~,aicL] = min(aic);
[~,bicL] = min(bic);
[~,hqcL] = min(hqc);

%Draw table if required
if strcmp(Table,'tableON')
display('-----------------------------------------------------------')
for i = 1:pmax
    AICstr='    bic=';
    BICstr='    hqc=';
    HQCstr=' ';
    if i==aicL
        AICstr='*   bic=';
    end
    if i==bicL
        BICstr='*   hqc=';
    end
    if i==hqcL
        HQCstr='*';    
    end
    display(['Lag ' num2str(i) ' aic=' num2str(aic(i)) AICstr num2str(bic(i)) BICstr num2str(hqc(i)) HQCstr ])
end
display('-----------------------------------------------------------')
end
end

% ANCILLARY FUNCTION
function [Y,X,Y_start,Y_hat,PI,e,cons]=varestimate(Data,p,constant)
% varestimate makes use of the SUR representation.
% INPUTS: 
% Data - t*n dataset at hand
% p - number of lags 
% constant - ==1 (==0) includes (excludes) the constant in the VAR
%
% OUTPUTS:
% Y - T*n matrix of the SUR representation
% X - T*(n*p) matrix of the SUR representation
% Y_start - p*n first p observations for each variable
% Y_hat - T*n fitted values
% PI - n*(n*p) matrix of coefficients [A1 A2 ... Ap]
% e - T*n estimated residuals
% cons - n*1 constant in the model (set =0 if constant==0)

% STEP 1: Generate Matrixies for SUR Representation
Y=(Data(p+1:end,:));          % lose first p lags 
X = lagmatrix(Data,1:p);      % generate lagged values 
X(1:p,:)=[];                  % Lose one osb for each lag
 
if constant==1;
  X=[ones(size(X,1),1) X];    % add the constant if needed
end
  
% Isolate the p inital elements of Data
Y_start=(Data(1:p,:)); 

% STEP 2: estimate the model 
MLE=((X'*X)\(X'*Y))'; 
if constant ==1;
cons = MLE(:,1);                % Isolate the constant
PI = MLE(:,2:end);              % Isolate the coefficients 
elseif constant ==0;
cons=zeros(size(MLE,1),1);      % set the constant = 0
PI = MLE;
end

% STEP 3: find fitted values and residuals:
Y_hat=X*MLE';
e = Y-Y_hat;
end

