function [BL]=computeirfs(BL,transf)
% This function compute irfs for factor model
% Inputs:
% - BL: impulse responses
% - transf: vector with index for transformation:
%          * 0: var in levels
%          * 1: var in first differences
%          * 2: var in log first differences
%          * 3: var in second difference
% Output: BL will be cumulated for first difference, cumulated and
% expressed in % points for log first diff and cumulated twice for second
% differences

ddf=find(transf==3);
df=find(transf==1);
dlf=find(transf==2);
BL(ddf,:,:)=cumsum(cumsum(BL(ddf,:,:),3),3);
BL(df,:,:)=cumsum(BL(df,:,:),3);
BL(dlf,:,:)=cumsum(BL(dlf,:,:),3)*100;

end