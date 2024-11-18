function y = objective(X,e_num,m,n,datanumber)
[e_num,m,n,band,o2,o2_3d,~,~,~] = todata(datanumber);
EM=[];
for a=1:e_num
   EM = [EM,o2(:,X(a))];
end
%% UCLS
s=hyperUcls(o2,EM);
s_new=max(0,s);
rec1=EM*s_new;
y(1)= rmse(o2,rec1,m,n,band);
%% FCLS
[~,rec2] = abundance(o2,m,n,EM);
y(2)= rmse(o2,rec2,m,n,band);
end