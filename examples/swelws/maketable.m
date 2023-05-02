gg = size(L);
Table(:,1)=L(:,1);
for jj=1:gg(2)-1;
    Table(:,jj*2)=L(:,jj+1);
    Table(:,jj*2+1)=0.01*round(100*(log2(L(:,jj)./L(:,jj+1))));
end