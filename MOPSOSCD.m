function [ps,pf]=MOPSOSCD(func_name,VRmin,VRmax,n_obj,Particle_Number,Max_Gen,e_num,m,n,datanumber,prb)
n_var=size(VRmin,2);
Max_FES=Max_Gen*Particle_Number;
n_PBA=5;
n_NBA=3*n_PBA;
cc=[2.05 2.05];                    %Acceleration constants
iwt=0.7298;                        %Inertia weight
%% Initialize particles' positions and velocities
    mv=0.5*(VRmax-VRmin);
    VRmin=repmat(VRmin,Particle_Number,1);
    VRmax=repmat(VRmax,Particle_Number,1);
    Vmin=repmat(-mv,Particle_Number,1);
    Vmax=-Vmin;
%% Initialization
vel = x_updata(Particle_Number,e_num,m,n);
pos = x_updata(Particle_Number,e_num,m,n);

fitness=zeros(Particle_Number,n_obj);
for i=1:Particle_Number
    fitness(i,:)=feval(func_name,pos(i,:),e_num,m,n,datanumber);
end
fitcount=Particle_Number;
particle=[pos,fitness];
row_of_cell=ones(1,Particle_Number); 
col_of_cell=size(particle,2);
PBA=mat2cell(particle,row_of_cell,col_of_cell);
NBA=PBA; 
for i=1:Max_Gen
    C=mean(pos(:,1:n_var));
    for j=1:Particle_Number
        if j==1
            tempNBA=PBA{Particle_Number,:};
            tempNBA=[tempNBA;PBA{1,:}];
            tempNBA=[tempNBA;PBA{2,:}];
        elseif j==Particle_Number
            tempNBA=PBA{Particle_Number-1,:};
            tempNBA=[tempNBA;PBA{Particle_Number,:}];
            tempNBA=[tempNBA;PBA{1,:}];
        else
            tempNBA=PBA{j-1,:};
            tempNBA=[tempNBA;PBA{j,:}];
            tempNBA=[tempNBA;PBA{j+1,:}];
        end
        NBA_j=NBA{j,1};
        tempNBA=[tempNBA;NBA_j];
        tempNBA=non_domination_scd_sort(tempNBA(:,1:n_var+n_obj), n_obj, n_var,datanumber);
        if size(tempNBA,1)>n_NBA
            NBA{j,1}=tempNBA(1:n_NBA,:);
        else
            NBA{j,1}=tempNBA;
        end
    end
    for k=1:Particle_Number

        updata_num = ceil(prb*Particle_Number);
        diff_order = randperm(Particle_Number);
        change_num = diff_order(:,(1:updata_num));
        PBA_k=PBA{k,1};
        pbest=PBA_k(1,:);
        NBA_k=NBA{k,:};
        nbest=NBA_k(1,:);
        % Update velocities
        vel(k,:)=iwt.*vel(k,:)+cc(1).*rand(1,n_var).*(pbest(1,1:n_var)-pos(k,:))+cc(2).*rand(1,n_var).*(nbest(1,1:n_var)-pos(k,:));
        % Make sure that velocities are in the setting bounds.
        vel(k,:)=(vel(k,:)>mv).*mv+(vel(k,:)<=mv).*vel(k,:);
        vel(k,:)=(vel(k,:)<(-mv)).*(-mv)+(vel(k,:)>=(-mv)).*vel(k,:);
        
        change_tf = ismember(k,change_num);
        if change_tf==0
            % Update positions
            pos(k,:)= ceil(pos(k,:)+vel(k,:));
        else
            pos(k,:) = pos(k,:);
        end
        for p=1:e_num
            if pos(k,p)<1 || pos(k,p)>m*n
                pos(k,p)=randperm(m*n,1);
            end
        end
        fitness(k,:)=feval(func_name,pos(k,1:n_var),e_num,m,n,datanumber);
        fitcount=fitcount+1;
        particle(k,1:n_var+n_obj)=[pos(k,:),fitness(k,:)];
        %% PBA
        PBA_k=[PBA_k(:,1:n_var+n_obj);particle(k,:)];
        PBA_k = non_domination_scd_sort(PBA_k(:,1:n_var+n_obj), n_obj, n_var,datanumber);
        if size(PBA_k,1)>n_PBA
            PBA{k,1}=PBA_k(1:n_PBA,:);
        else
            PBA{k,1}=PBA_k;
        end
    end
    
    if fitcount>Max_FES
        break;
    end
    global_ucls(i)= PBA_k(1,(e_num+1));
    global_fcls(i)= PBA_k(1,(e_num+2));
end

%% Output ps and pf
tempEXA=cell2mat(NBA);
tempEXA=non_domination_scd_sort(tempEXA(:,1:n_var+n_obj), n_obj, n_var,datanumber);
if size(tempEXA,1)>Particle_Number
    EXA=tempEXA(1:Particle_Number,:);
else
    EXA=tempEXA;
end
tempindex=find(EXA(:,n_var+n_obj+1)==1);
ps=EXA(tempindex,1:n_var);
pf=EXA(tempindex,n_var+1:n_var+n_obj);
end



function f = non_domination_scd_sort(x, n_obj, n_var,datanumber)

    [N_particle, ~] = size(x);
    front = 1;
% There is nothing to this assignment, used only to manipulate easily in MATLAB.
    F(front).f = []; 
    individual = [];
%% Non-Dominated sort. 
    for i = 1 : N_particle
        individual(i).n = 0; 
        individual(i).p = [];
        for j = 1 : N_particle
            dom_less = 0;
            dom_equal = 0;
            dom_more = 0;
            for k = 1 : n_obj
                if (x(i,n_var + k) < x(j,n_var + k))
                    dom_less = dom_less + 1;
                elseif (x(i,n_var + k) == x(j,n_var + k))  
                    dom_equal = dom_equal + 1;
                else
                    dom_more = dom_more + 1;
                end
            end
            if dom_less == 0 && dom_equal ~= n_obj  
                individual(i).n = individual(i).n + 1;
            elseif dom_more == 0 && dom_equal ~= n_obj  
                individual(i).p = [individual(i).p j];
            end
        end   
        if individual(i).n == 0  
            x(i,n_obj + n_var + 1) = 1;
            F(front).f = [F(front).f i];
        end
    end

    while ~isempty(F(front).f)
       Q = [];
       for i = 1 : length(F(front).f)
           if ~isempty(individual(F(front).f(i)).p)
                for j = 1 : length(individual(F(front).f(i)).p)
                    individual(individual(F(front).f(i)).p(j)).n = ...
                        individual(individual(F(front).f(i)).p(j)).n - 1;
                    if individual(individual(F(front).f(i)).p(j)).n == 0
                        x(individual(F(front).f(i)).p(j),n_obj + n_var + 1) = ...
                            front + 1;
                        Q = [Q individual(F(front).f(i)).p(j)];
                    end
               end
           end
       end
       front =  front + 1;
       F(front).f = Q;
    end
    [~,index_of_fronts] = sort(x(:,n_obj + n_var + 1));
    for i = 1 : length(index_of_fronts)
        sorted_based_on_front(i,:) = x(index_of_fronts(i),:);
    end
    current_index = 0;

%% SCD. Special Crowding Distance

    for front = 1 : (length(F) - 1)
  
        crowd_dist_obj = 0;
        y = [];
        previous_index = current_index + 1;
        for i = 1 : length(F(front).f)
            y(i,:) = sorted_based_on_front(current_index + i,:); 
        end
        current_index = current_index + i;
        sorted_based_on_objective = [];
        for i = 1 : n_obj+n_var
            [sorted_based_on_objective, index_of_objectives] = sort(y(:,i));
            sorted_based_on_objective = [];
            for j = 1 : length(index_of_objectives)
                sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);
            end
            f_max = sorted_based_on_objective(length(index_of_objectives), i);
            f_min = sorted_based_on_objective(1,i);

            if length(index_of_objectives)==1
                y(index_of_objectives(1),n_obj + n_var + 1 + i) = 1;
            elseif i>n_var
                y(index_of_objectives(1),n_obj + n_var + 1 + i) = 1;
                y(index_of_objectives(length(index_of_objectives)),n_obj + n_var + 1 + i)=0;
            else
                 y(index_of_objectives(length(index_of_objectives)),n_obj + n_var + 1 + i)...
                    = 2*(sorted_based_on_objective(length(index_of_objectives), i)-...
                sorted_based_on_objective(length(index_of_objectives) -1, i))/(f_max - f_min+1);
                 y(index_of_objectives(1),n_obj + n_var + 1 + i)=2*(sorted_based_on_objective(2, i)-...
                sorted_based_on_objective(1, i))/(f_max - f_min+1);
            end
             for j = 2 : length(index_of_objectives) - 1
                next_obj  = sorted_based_on_objective(j + 1, i);
                previous_obj  = sorted_based_on_objective(j - 1,i);
                if (f_max - f_min == 0)
                    y(index_of_objectives(j),n_obj + n_var + 1 + i) = 1;
                else
                    y(index_of_objectives(j),n_obj + n_var + 1 + i) = ...
                         (next_obj - previous_obj)/(f_max - f_min);
                end
             end
        end
    %% Calculate distance in decision space
        crowd_dist_var = [];
        crowd_dist_var(:,1) = zeros(length(F(front).f),1);
        for i = 1 : n_var
            crowd_dist_var(:,1) = crowd_dist_var(:,1) + y(:,n_obj + n_var + 1 + i);
        end
        crowd_dist_var=crowd_dist_var./n_var;
        avg_crowd_dist_var=mean(crowd_dist_var);
    %% Calculate distance in objective space
        crowd_dist_obj = [];
        crowd_dist_obj(:,1) = zeros(length(F(front).f),1);
        for i = 1 : n_obj
            crowd_dist_obj(:,1) = crowd_dist_obj(:,1) + y(:,n_obj + n_var + 1+n_var + i);
        end
        crowd_dist_obj=crowd_dist_obj./n_obj;
        avg_crowd_dist_obj=mean(crowd_dist_obj);
    %% Calculate special crowding distance
        special_crowd_dist=zeros(length(F(front).f),1);
        for i = 1 : length(F(front).f)
            if crowd_dist_obj(i)>avg_crowd_dist_obj||crowd_dist_var(i)>avg_crowd_dist_var
                special_crowd_dist(i)=max(crowd_dist_obj(i),crowd_dist_var(i));
            else
                special_crowd_dist(i)=min(crowd_dist_obj(i),crowd_dist_var(i));
            end
        end
        y(:,n_obj + n_var + 2) = special_crowd_dist;
        y(:,n_obj + n_var + 3) = crowd_dist_var;
        y(:,n_obj + n_var + 4) = crowd_dist_obj;
        [~,index_sorted_based_crowddist]=sort(special_crowd_dist,'descend');
        y=y(index_sorted_based_crowddist,:);
        y = y(:,1 : n_obj + n_var+4 );
        z(previous_index:current_index,:) = y;
    end
    
f = z();
end