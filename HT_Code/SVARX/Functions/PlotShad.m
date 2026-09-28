function PlotShad(x,tim,color,lin)
if nargin==2
    color=1;
    lin=1;
    col2 = [.7 .7 .7];
elseif nargin>2
    if color==1 %light gray
        col2 = [.7 .7 .7];
    elseif color==2 %dark gray
        col2 = [.5 .5 .5];
    elseif color==3 %red
        col2 = [1 0 0];
    end
end

if nargin==1
    col2 = [.7 .7 .7];
    n=size(x,1);
    ll=10;
    TEMP = x(:,2);
    te4 = x(:,1);
    te5 = x(:,3);
    b=fill([1:n n:-1:1],[te4' flipud(te5)'],col2);
    set(b,'EdgeAlpha',0);
    hold on
    a = plot(1:n,[squeeze(TEMP) ]);
    set(a(1),'Color','k','LineWidth',1.5);
    z = plot(1:n,zeros(n,1));
    set(z(1),'Color','r','LineWidth',1);
    hold off
else
    n=size(x,1);
    ll=10;
    TEMP = x(:,2);
    te4 = x(:,1);
    te5 = x(:,3);
    b=fill([tim fliplr(tim)],[te4' flipud(te5)'],col2);
    set(b,'EdgeAlpha',0);
    hold on
    a = plot(tim,[squeeze(TEMP) ]);
    if lin==1
        set(a(1),'Color','k','LineStyle','-','LineWidth',1.5);
    elseif lin==2
        set(a(1),'Color','k','LineStyle','--','LineWidth',1.5);
    elseif lin==3
        set(a(1),'Color','k','LineStyle',':','LineWidth',1.5);
    end
    z = plot(tim,zeros(n,1));
    set(z(1),'Color','r','LineWidth',1);
    hold off
end
