function result = main(M,K,C,F,u0,ut0,t_final,dt,N_snap,nRandomVarsUsed,N,GS_ord,Kmodes,nint,remain,remove,dfix,Method_TM2,output)
fprintf('Start generating results...\n');

T = 0:dt:t_final;
T_snap=0:(t_final/N_snap):t_final; %snapshot time points
szt = length(T);

%load Legendre PC 
load('PC_Legendre.mat');
PC = PC_Legendre{nRandomVarsUsed};

%% MCS
fprintf('-------------------------------------------------------------\n');
fprintf('Monte Carlo simulation running...\n');
fprintf('-------------------------------------------------------------\n');
xisample=rand(nRandomVarsUsed,N); %Uniform sampling in (0,1)
xisample=2*xisample-1; %Uniform sampling in (-1,1)
tic;
[mean_MC,std_MC] = MCS(M,K,C,F,u0,ut0,szt,dt,Method_TM2,nRandomVarsUsed,xisample,remain,remove,dfix,N,output);
t_MC = toc;
fprintf('timeMC=%f\n', t_MC);
result.mean_MCS = mean_MC;
result.std_MCS = std_MC;
result.t_MCS = t_MC;

%% GS
fprintf('-------------------------------------------------------------\n');
fprintf('GS scheme running...\n');
fprintf('-------------------------------------------------------------\n');
for i = 1:GS_ord
    PC_i = PC{i};
    tic;
    show = 1;
    [mean_GS_all{i},var_GS_all{i}, ~, ~, ~] = GS(PC_i, M, K, C, F, ...
        nint, nRandomVarsUsed, T_snap, dt, u0,ut0,Method_TM2,remain,remove,dfix,show);
    mean_GS(i,:) = mean_GS_all{i}(output.id,:);
    std_GS(i,:) = var_GS_all{i}(output.id,:).^0.5;
    t_GS(i)=toc;
    string =['\norder = ',num2str(i), ' timeGS = ', num2str(t_GS(i)), '\n'];
    fprintf(string);
    fprintf('\n');
end
result.mean_GS = mean_GS;
result.std_GS = std_GS;
result.t_GS = t_GS;

%% GSD
fprintf('-------------------------------------------------------------\n');
fprintf('GSD scheme running...\n');
fprintf('------------------------------------------------------------\n');
for i = 1:length(Kmodes)
    fprintf('\n Kmodes=%i', Kmodes(i));
    tic;
    result_GSD{i} = GSD(PC{2}, M, K, C, F, Kmodes(i), ...
            nRandomVarsUsed, nint,dt,Method_TM2, T_snap,u0,ut0,remain,remove,dfix);
    t_GSD(i)=toc;
    iter = result_GSD{i}.iter;
    mean_GSD(i,:) = result_GSD{i}.mean{iter}(output.id,:);
    std_GSD(i,:) = result_GSD{i}.var{iter}(output.id,:).^0.5;
end
result.mean_GSD = mean_GSD;
result.std_GSD = std_GSD;
result.t_GSD = t_GSD;

%% AAPG
fprintf('-------------------------------------------------------------\n');
fprintf('AAPG scheme running...\n');
fprintf('-------------------------------------------------------------\n');
AAPG_GS_ord = 2; 
[mean_AAPG, std_AAPG,t_AAPG]= AAPG(PC{AAPG_GS_ord}, M, K, C, F, nint, ...
    nRandomVarsUsed, AAPG_GS_ord, PC_Legendre{1}{AAPG_GS_ord},PC_Legendre{2}{AAPG_GS_ord}, szt, T_snap, dt, u0, ut0, Method_TM2, remain, remove, dfix, output);
result.mean_AAPG = mean_AAPG;
result.std_AAPG = std_AAPG;
result.t_AAPG = t_AAPG;
