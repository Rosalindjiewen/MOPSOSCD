function O=FCLS(Data,A)
[p,m,n]=size(Data);
[p,classnum]=size(A);
O=zeros(m,n,classnum);
for i=1:m
    for j=1:n
        if all(Data(:,i,j)==0)
            Data(1,i,j)=0.0001;
        end
        abundance_test= FCLS_single_pixel(A,Data(:,i,j));
        for k=1:classnum
            O(i,j,k) = abundance_test(k);
        end
    end
end
