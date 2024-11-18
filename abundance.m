function [s_new,rec]=abundance(V,m,n,E)
X_3d = hyperConvert3d(V, m, n);
[r_i,r_j,r_k]=size(X_3d);
for i=1:r_i
    for j=1:r_j
        for k=1:r_k
            X_3d_new(k,i,j)=X_3d(i,j,k);
        end
    end
end  
s=FCLS(X_3d_new,E);
s_new=hyperConvert2d(s);
rec=E*s_new;
end