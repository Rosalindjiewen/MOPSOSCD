clear all
clc

global fname
%% parameters
N_function=1; 
%data
datanumber = 4; 
[e_num,m,n,band,o2,o2_3d,m_turth,filename,RGBband] = todata(datanumber);
save_dir = 'ans/';   %save_dir
popsize = 3;
Max_Gen = 9;
prb = 0.2;
for i=1:N_function
    switch i
        case 1
            fname='objective';   
            n_obj=2;
            n_var=e_num;
            xl=[];
            xu=[];
            for k=1:n_var
                xl=[xl,0];
                xu=[xu,m*n];
            end
    end
   fprintf('Running test function: %s \n', fname);
   %% Search the PSs using MOPSOSCD
    [ps,pf]= MOPSOSCD(fname,xl,xu,n_obj,popsize,Max_Gen,e_num,m,n,datanumber,prb);
        
end
%% save_ps & pf
n_ps = strcat(save_dir,'ps.mat');
save(n_ps,'ps');
n_pf = strcat(save_dir,'pf.mat');
save(n_pf,'pf');
