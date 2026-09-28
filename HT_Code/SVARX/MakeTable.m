function Table = MakeTable(Table,ndec,rowtitles, columntitles)
[I,J]=size(Table);
if nargin==3
    
    for i=1:J
    columntitles{i}='';
    end
end
if nargin==2
    for i=1:I
    rowtitles{i}='';
    end
    for i=1:J
    columntitles{i}='';
    end
end
Table=round(Table*10^ndec)/10^ndec;


for i=1:I
    fprintf(rowtitles{i})
for j=1:J
    if ndec==2
    fprintf('  & %8.2f ', Table(i,j))
    elseif ndec==1
     fprintf(' & %8.1f ', Table(i,j)) 
     elseif ndec==3
     fprintf(' & %8.3f ', Table(i,j)) 
     elseif ndec==0
     fprintf('  & %8.0f ', Table(i,j))
    else
     fprintf(' & %8.4f ', Table(i,j)) 
    end
end
    fprintf('\\\\ \n')
    
end
