function [abundance]= FCLS_single_pixel(A,r1)
% A: endmember matrix, r1: mixture
% tic
%Dan's: optimal
numloop=size(A,2);  % Class Number
e=1/10/max(r1);  % e is delta in the algorithm
eA=e*A;
E=[ones(1,numloop);eA];
EtE=E'*E;
[m,n] = size(EtE);
One=ones(m,1);
%if rank(EtE)==size(EtE,1)
    iEtE=pinv(EtE);
%else
    %iEtE=invertbysvd(EtE);
%end
iEtEOne=iEtE*One;  % s
sumiEtEOne=sum(iEtEOne); % 1's
weights=diag(iEtE);

c=0;
sample=r1;
er=e*sample;        
f=[1;er];
Etf=E'*f;

tol=1e-7;
%fcls1a

%%%% THIS IS lamdiv2
ls=iEtE*Etf;
lamdiv2=(1-(ls'*One))/sumiEtEOne;        % ??
x2=ls-lamdiv2*iEtEOne;
x2old=x2;
iter_on_whole=0;
if (any(x2<-tol))&&(iter_on_whole<=10*m)
    iter_on_whole=iter_on_whole+1;
    Z=zeros(m,1);
    iter=0;
    while(any(x2<-tol) && iter <(m))
        Z(x2<-tol)=1;
        zz=find(Z);
        x2=x2old;              % Reset x2
        L=iEtE(zz,zz);
        ab=size(zz);
        lastrow=ab(1)+1;   
        lastcol=lastrow;
        L(lastrow,1:ab(1))=(iEtE(:,zz)'*One)';
        L(1:ab(1),lastcol)=iEtEOne(zz);
        L(lastrow,lastcol)=sumiEtEOne;
        xerow=x2(zz);
        xerow(lastrow,1)=0;
        lagra=lsqminnorm(L,xerow);
        while (any(lagra(1:ab(1))>0))
            maxneg=weights(zz).*lagra(1:ab(1));
            [yz,iz]=max(maxneg);
            Z(zz(iz))=0;
            zz=find(Z);
            L=iEtE(zz,zz);
            ab=size(zz);
            lastrow=ab(1)+1;
            lastcol=lastrow;
            L(lastrow,1:ab(1))=(iEtE(:,zz)'*One)';
            L(1:ab(1),lastcol)=iEtEOne(zz);
            L(lastrow,lastcol)=sumiEtEOne;
            xerow=x2(zz);    
            xerow(lastrow,1)=0;
            lagra=lsqminnorm(L,xerow);
        end
        %problem with lamscls zz may be null
        if ~isempty(zz)
            x2=x2-iEtE(:,zz)*lagra(1:ab(1))-lagra(lastrow)*iEtEOne;
        end
        iter=iter+1;
    end
end
abundance=x2;


