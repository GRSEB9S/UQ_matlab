function [phi_coeffs,dphi_coeffs,ddphi_coeffs] = GSD_secondsubproblemP2(PC_lambda, massmat, stiff, damp, source, nRandomVarsUsed, nint, T_snap, dt, initial_dis, initial_vel, Method_TM2,remain,remove,dfix)
show = 1;
[~, ~, phi_coeffs_mat,dphi_coeffs_mat, ddphi_coeffs_mat] = GS(PC_lambda, massmat, stiff, damp, source, nint, nRandomVarsUsed, T_snap, dt, initial_dis,initial_vel,Method_TM2,remain,remove,dfix,show);
for i=1:PC_lambda.Size
   phi_coeffs{i}=phi_coeffs_mat((i-1)*nint+1:i*nint,:);
   dphi_coeffs{i}=dphi_coeffs_mat((i-1)*nint+1:i*nint,:);
   ddphi_coeffs{i}=ddphi_coeffs_mat((i-1)*nint+1:i*nint,:);
end