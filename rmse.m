function y = rmse(r,rec,m,n,row)
sub_r = r-rec;
for t = 1:m*n
    k1(:,t) = norm(sub_r(:,t),2);
    k2(:,t) = sqrt((1/row)*((k1(:,t))^2));
end
y = (1/(m*n)*sum(k2));
end