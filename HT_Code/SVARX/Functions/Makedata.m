% Import and name data:

%[data,header]=xlsread('GMS_data.xls','Monthly','B2:O589');

load 'data'
load 'RmRf'

ip=log(data(:,1))*100;
cpi=log(data(:,2))*100;
ebp=data(:,3);
ffr=data(:,4);
sp500=log(data(:,5))*100;
pce=log(data(:,6)./(data(:,7)./100))*100;
pcepilfe=log(data(:,7)).*100;
finstress=log(data(:,8))*100;
finuncertainty=data(:,9);
macrouncertainty=data(:,14);
unrate=data(:,10);
NFCI=data(:,11);
GS10=data(:,12);
GS10state=data(:,13);
%RmRf=data(:,15);