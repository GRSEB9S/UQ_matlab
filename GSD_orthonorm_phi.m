%
% Orthogonalization of functions \Phi(t)
%

function [phi_coeffs_out] = GSD_orthonorm_phi(dt_euler, Nt, t_final, Kmodes, phi_coeffs_in)

%Gram-Schmidt procedure:
ps11sq=GSD_pscal_phi(phi_coeffs_in{1}, phi_coeffs_in{1}, dt_euler, Nt, t_final);
phi_coeffs_out{1}=phi_coeffs_in{1}/sqrt(ps11sq);

for j=2:Kmodes
   tmp=phi_coeffs_in{j};
   for k=1:j-1
       tmp=tmp-GSD_pscal_phi(phi_coeffs_in{j}, phi_coeffs_out{k}, dt_euler, Nt, t_final)*phi_coeffs_out{k};
   end
   %normalization:
   psjjsq=GSD_pscal_phi(tmp, tmp, dt_euler, Nt, t_final);
   phi_coeffs_out{j}=tmp/sqrt(psjjsq);
end
