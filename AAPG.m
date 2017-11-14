function [mean_AAPG, std_AAPG,t_AAPG] = AAPG(PCLegendre, M, K, C, F, nint, ...
    nRandomVarsUsed, orderPolSol, PC_1, PC_2, Nt, T_snap, dt_euler, init_dis, init_vel, Method_TM2, remain, remove, dfix, output)
%No. of sub-problems in first and second order AAPG
PA(1) = 1+nRandomVarsUsed; %first-order
PA(2) = PA(1)+nRandomVarsUsed*(nRandomVarsUsed-1)/2; %second order
    
%compute zero-order component function
[u0,~,~] = Second_solver(Method_TM2,init_dis,init_vel,dt_euler,Nt,K{1},C{1},M{1},F{1},0);

%compute first-order component functions
tic;
fprintf('Computing first-order components...\n');
uprime{1}=first_order(PC_1, M, K, C, F,nRandomVarsUsed, ...
    T_snap, dt_euler, nint, init_dis, init_vel, Method_TM2, remain, remove, dfix);
fprintf('Assembling results...\n');
%assemble results
for j = 1:nint
    coeff{j} = ones(PA(1),Nt);
end
Solu{1} = assembleANOVA(nRandomVarsUsed,orderPolSol,u0,uprime,1,nint/2,coeff);
t_AAPG(1)=toc;
string =['time_AAPG1=',num2str(t_AAPG(1)),'\n'];
fprintf(string);

%compute second-order component functions
tic;
fprintf('Computing second-order components...\n');
uprime{2}=second_order(PC_2, M, K, C, F, ...
    nRandomVarsUsed, T_snap, dt_euler, nint, init_dis, init_vel, Method_TM2, remain, remove, dfix);
fprintf('Assembling results...\n');
for j = 1:nint
    coeff{j} = ones(PA(2),Nt);
end
Solu{2} = assembleANOVA(nRandomVarsUsed,orderPolSol,u0,uprime,2,nint/2,coeff);
t_AAPG(2)=toc;
string =['time_AAPG2=',num2str(t_AAPG(2)),'\n'];
fprintf(string);

for i = 1:2
    %mean and variance of the solution
    m_dis_AAPG{i}=Solu{i}{1};
    index_retain=ones(size(Solu{i},1),1);
    v_dis_AAPG{i}=get_variance(Solu{i}, PCLegendre, nint, Nt, index_retain);

    %apply the boundary conditions
    mean_AAPG_all{i} = fixconstrain(m_dis_AAPG{i},'Dis',remain,remove,dfix);
    var_AAPG_all{i} = fixconstrain(v_dis_AAPG{i},'Var',remain,remove,dfix);
    
    %fetch specific dof
    mean_AAPG(i,:) = mean_AAPG_all{i}(output.id,:);
    std_AAPG(i,:) = var_AAPG_all{i}(output.id,:).^0.5;
end
end