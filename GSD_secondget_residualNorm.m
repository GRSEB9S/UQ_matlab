function [normResid] = GSD_secondget_residualNorm(PC, phi_coeffs, dphi_coeffs, ddphi_coeffs, lambda_coeffs, massmat, stiff, damp, source, Kmodes, nRandomVarsUsed, nint, dt, T_snap)
t_final = T_snap(length(T_snap));

sztime=size(phi_coeffs{1},2);
PsiSqNorm=PC.PsiSqNorm;

for m=1:nRandomVarsUsed+1
   for k=1:Kmodes
      Mass_ddphi{m}{k}=massmat{m}*ddphi_coeffs{k};
      Stiff_phi{m}{k}=stiff{m}*phi_coeffs{k};
      Damp_dphi{m}{k}=damp{m}*dphi_coeffs{k};
   end
end

for j=1:PC.Size
   Resid{j}=zeros(nint,sztime);
   tmp=zeros(nint,sztime);
   for k=1:Kmodes
     for p=1:PC.Size
        for m=1:nRandomVarsUsed+1
             tmp=tmp+PC.DD{m}(p,j)*lambda_coeffs{k}(p)*(Mass_ddphi{m}{k}+Stiff_phi{m}{k}+Damp_dphi{m}{k});
        end
     end
   end
   Resid{j}=Resid{j}+tmp;
   for m = 1:nRandomVarsUsed+1
       Resid{j}=Resid{j}-source{m}*PC.Psilambda(m,j);
   end
end

% for k=1:PC.Size
%   MassResid{k}=massmat*Resid{k};
% end
normResid=0.0;
for k=1:PC.Size
   normResid=normResid+PsiSqNorm(k)*dt*sum(diag((Resid{k})'*Resid{k}));
end