 function y=x_updata(NP,e_num,m,n)
  y=[];
  for i=1:NP
      a=randperm(m*n,e_num);
      y=[y;a];
  end
  end