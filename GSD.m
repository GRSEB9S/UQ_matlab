function result = GSD(PC, massmat, stiff, damp, source, Kmodes, ...
  nRandomVarsUsed, nint,dt,Method_TM2,T_snap,initial_dis,initial_vel,remain,remove,dfix)

Nt = T_snap(length(T_snap))/dt+1;
t_final = T_snap(length(T_snap));

time0 = cputime;

%---stochastic coefficients;
lambda_coeffs=cell(Kmodes,1);
for i=1:Kmodes
  lambda_coeffs{i}=rand(PC.Size,1); %random initialization here
end
%---reorthogonalization:
lambda_coeffs=GSD_orthonorm_lbd(PC, Kmodes, lambda_coeffs);


%--------------------------------------MAIN LOOP-------------------------------------------------%

iter_max=15;
iter=1;
tol=1e-8; %tolerance for stopping criterion
err_iter=1.0;
maxerr_var_2iters=1.0;
diff_normResid=1.0;

while ( (iter<=iter_max) & (diff_normResid>tol) )

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPROBLEM P2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %---compute stochastic elements related to the lambda_coeffs:
     PC_lambda=GSD_get_stoch_lambda(PC, lambda_coeffs, Kmodes, nRandomVarsUsed);
     %---compute space-time coefficients:
     [phi_coeffs,dphi_coeffs, ddphi_coeffs]=GSD_secondsubproblemP2(PC_lambda, massmat, stiff, damp, source, nRandomVarsUsed, nint, T_snap, dt, initial_dis, initial_vel, Method_TM2 ,remain,remove,dfix);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     phi_coeffs_stock{iter}=phi_coeffs;
     dphi_coeffs_stock{iter}=dphi_coeffs;
     lambda_coeffs_stock{iter}=lambda_coeffs;

     %COMPUTE STATISTICS OF GSD:
     [muGSD, varGSD]=GSD_get_statistics(phi_coeffs, lambda_coeffs, PC, Kmodes, nint, Nt);
     varGSD_stock{iter}=varGSD;
     meanGSD_stock{iter}=muGSD;

     %COMPUTE RESIDUAL NORM:
     normResid=GSD_secondget_residualNorm(PC, phi_coeffs, dphi_coeffs, ddphi_coeffs, lambda_coeffs, massmat, stiff, damp, source, Kmodes, nRandomVarsUsed, nint, dt, T_snap);
     normResid_stock{iter}=normResid;
 
     %---reorthogonalization:
     phi_coeffs=GSD_orthonorm_phi(dt, Nt, t_final, Kmodes, phi_coeffs);

     %COMPUTE ERROR BETWEEN TWO SUCCESSIVE GSD EXPANSIONS:
     if (iter >=2)
        err_iter=GSD_get_err_iter(phi_coeffs_stock{iter}, lambda_coeffs_stock{iter}, ...
                              phi_coeffs_stock{iter-1}, lambda_coeffs_stock{iter-1}, PC, Kmodes, nint, Nt, dt);
     end

     if (iter >=2)
        diff_normResid(iter-1)=abs(normResid_stock{iter}-normResid_stock{iter-1});
        string=['diff_normResid=', num2str(diff_normResid(iter-1)),'\n'];
        fprintf(string);
     end

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPROBLEM P1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %---compute stochastic coefficients;
     [lambda_coeffs] = GSD_secondsubproblemP1(PC, massmat, stiff, damp, source, Kmodes, phi_coeffs, dphi_coeffs, ddphi_coeffs, nRandomVarsUsed, dt);
     %---reorthogonalization:
     lambda_coeffs=GSD_orthonorm_lbd(PC, Kmodes, lambda_coeffs);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     iter=iter+1;

end
time = cputime - time0;
iter_out=iter-1;
%----------------------------------END OF MAIN LOOP-------------------------------------------------%

for i = 1:iter_out
    meanGSD_stock{i} = fixconstrain(meanGSD_stock{i},'Dis',remain,remove,dfix);
    varGSD_stock{i} = fixconstrain(varGSD_stock{i},'Var',remain,remove,dfix);
end
result.mean = meanGSD_stock;
result.var = varGSD_stock;
result.iter = iter_out;
result.resid = normResid_stock;
result.diffresid = diff_normResid;
result.time = time;
phi_coeffs=phi_coeffs_stock{iter-1};
dphi_coeffs=dphi_coeffs_stock{iter-1};
lambda_coeffs=lambda_coeffs_stock{iter-1};
beta_GSD = zeros(PC.Size*nint,Nt);
dbeta_GSD = zeros(PC.Size*nint,Nt);
for k=1:PC.Size
    tmp=zeros(nint,Nt);dtmp=zeros(nint,Nt);
    for j=1:Kmodes
        tmp=tmp+lambda_coeffs{j}(k)*phi_coeffs{j};
        dtmp=dtmp+lambda_coeffs{j}(k)*dphi_coeffs{j};
    end
    beta_GSD((k-1)*nint+1:k*nint,:)=tmp;
    dbeta_GSD((k-1)*nint+1:k*nint,:)=dtmp;
end
result.sol_dis_2GSD = beta_GSD;
result.sol_vel_2GSD = dbeta_GSD;
fprintf('\n');
fprintf('\n timeGSD=%f', time);
fprintf('\n iter_GSD=%i\n', iter_out);
fprintf('\n');
