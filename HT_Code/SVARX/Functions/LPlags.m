%%% Function to perform Local Projections a la Jordà (2005)
%%% Author: Nicolo' Maffei Faccioli

function [lp,bands,errors,vardecdenom]=LPlags(y,lags,controls,X,hor,c)

    p=hor;
    errors=NaN(length(y),p+1);
    lp=zeros(size(controls,2)+c+lags+size(X,2),p+1);
    bands=zeros(size(controls,2)+c+lags+size(X,2),p+1);
    vardecdenom=zeros(1,p+1);
    lagdep=lagmaker(y,lags); % lagged dependent variable
    
    for i=1:p+1
    
    if lags>=i
        dep=y(i+lags+1-i:end);
        regX=lagmatrix(X,i-1);
        regX=regX(i+lags+1-i:end,:);
        reg=[lagdep controls(i+lags+1-i:end,:) regX];
    else, dep=y(i+1-1:end);
          regX=lagmatrix(X,i-1);
          regX=regX(i+1-1:end,:);
          reg=[lagdep(i-lags:end,:) controls(i+1-1:end,:) regX];
    end
    
    [~,bands(:,i),lp(:,i)]=hac(reg,dep,'display','off');
    errors=dep-[ones(length(dep),1),reg]*lp(:,i);
    
    % LPA:
    % vardecdenom(:,i)=var((lp(c+lags+size(controls,2)+1,i).^2).*regX+errors);
    
    % LPB:
    
    if i==1
    vardecdenom(:,i)=var(errors);
    else 
        
        if lags>=i
            regX=lagmaker(X,i-1);
            vardecdenom(:,i)=var(errors-regX(lags:end,:)*(lp(c+lags+size(controls,2)+1,1:(i-1)))');
        else, vardecdenom(:,i)=var(errors-lagmaker(X,i-1)*(lp(c+lags+size(controls,2)+1,1:(i-1)))');
        
        end
    end
    
    end
    
end
