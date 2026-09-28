Shocknames={'Financial Shock Squared'};

magnitude1=1;

C_magnitude1=magnitude1.*C_wold(:,2:m:end);

HighC_magnitude1=HighC(:,2:m:end);
LowC_magnitude1=LowC(:,2:m:end);
HighC90_magnitude1=HighC90(:,2:m:end);
LowC90_magnitude1=LowC90(:,2:m:end);

magnitude2=2;

C_magnitude2=(magnitude2.^2).*C_wold(:,2:m:end);
%C_magnitude2=(magnitude2).*C_wold(:,2:m:end);% absolute value

magnitude3=4;

C_magnitude3=(magnitude3.^2).*C_wold(:,2:m:end);
%C_magnitude3=(magnitude3).*C_wold(:,2:m:end); % absolute value

irf_plot_VARX_magnitude(C_magnitude1,HighC90_magnitude1,LowC90_magnitude1,HighC_magnitude1,LowC_magnitude1,C_magnitude2,C_magnitude3,VARnames,Shocknames,colorBNDS,colorBNDS68,hor,n,m)
legend({'90% confidence bands','68% confidence bands','IRF','IRF 2 sd','IRF 4 sd'},'FontSize',18,'Orientation','Horizontal')