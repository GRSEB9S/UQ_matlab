function [lambda_coeffs] = GSD_secondsubproblemP1(PC, massmat, stiff, damp, source, Kmodes, phi_coeffs, dphi_coeffs, ddphi_coeffs, nRandomVarsUsed,dt)

K_phi=cell(nRandomVarsUsed+1,Kmodes);
M_ddphi=cell(nRandomVarsUsed+1,Kmodes);
C_dphi=cell(nRandomVarsUsed+1,Kmodes);
for m=1:nRandomVarsUsed+1
    for j=1:Kmodes
       K_phi{m}{j}=stiff{m}*phi_coeffs{j};
       M_ddphi{m}{j}=massmat{m}*ddphi_coeffs{j};
       C_dphi{m}{j}=damp{m}*dphi_coeffs{j};
    end
end

%---compute integrals in time:
A=zeros(Kmodes,Kmodes);B=zeros(Kmodes,Kmodes);C=zeros(Kmodes,Kmodes);
for i=1:Kmodes
    for j=1:Kmodes
        A(j,i)=dt*sum(diag(M_ddphi{1}{i}*phi_coeffs{j}'));
        B(j,i)=dt*sum(diag(K_phi{1}{i}*phi_coeffs{j}'));
        C(j,i)=dt*sum(diag(C_dphi{1}{i}*phi_coeffs{j}'));
    end
end

D=cell(nRandomVarsUsed,1);
E=cell(nRandomVarsUsed,1);
F=cell(nRandomVarsUsed,1);
for m=2:nRandomVarsUsed+1
    for i=1:Kmodes
      for j=1:Kmodes
         D{m-1}(j,i)=dt*sum(diag(M_ddphi{m}{i}*phi_coeffs{j}'));
         E{m-1}(j,i)=dt*sum(diag(K_phi{m}{i}*phi_coeffs{j}'));
         F{m-1}(j,i)=dt*sum(diag(C_dphi{m}{i}*phi_coeffs{j}'));
      end
    end
end

%---assemble rhs:
source_phi_coeffs=cell(nRandomVarsUsed+1,1);
for m = 1:nRandomVarsUsed+1
    source_phi_coeffs{m}=zeros(Kmodes,1);
    for i=1:Kmodes
        source_phi_coeffs{m}(i)=dt*sum(diag(source{m}*phi_coeffs{i}'));
    end
end

%solve the stochastic linear system using Prasanth's code:
AA=cell(nRandomVarsUsed+1);
AA{1}=A+B+C;
for m=2:nRandomVarsUsed+1
    AA{m}=D{m-1}+E{m-1}+F{m-1};
end
bb=cell(nRandomVarsUsed+1,1);
for m = 1:nRandomVarsUsed+1
    bb{m}=source_phi_coeffs{m};
end

[~, ~, xcell]=GhanemSpanos_linearSystem(AA, bb, PC.PsiSqNorm, PC.DD, PC.Size);
for i=1:Kmodes
   for p=1:PC.Size
      lambda_coeffs{i}(p,1)=xcell{p}(i);
   end
end