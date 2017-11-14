function [mu_dis, var_dis, soldis, solvel, solacc] = GS(PC, massmat, stiff, damp, source, nint, nRandomVarsUsed, T_snap, dt, initial_dis,initial_vel,Method_TM2,remain,remove,dfix,show)
t_final = T_snap(length(T_snap));

P=PC.Size;
PsiSqNorm=PC.PsiSqNorm;
Nt=t_final/dt+1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Kbig,Cbig,Mbig and Fbig        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kbig = sparse(P*nint,P*nint);
Mbig = sparse(P*nint,P*nint);
Cbig = sparse(P*nint,P*nint);

Fbig = sparse(nint*P,Nt);
for m = 1:nRandomVarsUsed+1
    Kbig = Kbig + kron(PC.DD{m},stiff{m});
    Mbig = Mbig + kron(PC.DD{m},massmat{m});
    Cbig = Cbig + kron(PC.DD{m},damp{m});
    for i=1:P
        Fbig((i-1)*nint+1:i*nint,:)=Fbig((i-1)*nint+1:i*nint,:)+PC.Psilambda(m,i)*source{m};
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      ODE resolution                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialization
dis0=zeros(P*nint,1);
vel0=zeros(P*nint,1);
for i = 1:P
    dis0((i-1)*nint+1:i*nint)=initial_dis*PC.MeanBasis(i)/PC.PsiSqNorm(i);
    vel0((i-1)*nint+1:i*nint)=initial_vel*PC.MeanBasis(i)/PC.PsiSqNorm(i);
end
dis0 = sparse(dis0);
vel0 = sparse(vel0);

[soldis,solvel,solacc] = Second_solver(Method_TM2,dis0,vel0,dt,Nt,Kbig,Cbig,Mbig,Fbig,show);

% Postprocessing:
mu_dis=soldis(1:nint,:); %mean of displacement
mu_vel=solvel(1:nint,:); %mean of velocity
var_dis=zeros(nint,Nt);var_vel=zeros(nint,Nt);
for i = 2:P
    dis_temp = soldis(nint*(i-1)+1:nint*i,:);
    var_dis = var_dis + dis_temp.^2*PsiSqNorm(i);
    vel_temp = solvel(nint*(i-1)+1:nint*i,:);
    var_vel = var_vel + vel_temp.^2*PsiSqNorm(i);
end
%output:
mu_dis = fixconstrain(mu_dis,'Dis',remain,remove,dfix);
var_dis = fixconstrain(var_dis,'Var',remain,remove,dfix);



