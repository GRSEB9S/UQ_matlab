%
% Inner product in L^2(0,T)
%

function [pscal] = GSD_pscal_phi(phi_i, phi_j, dt, Nt, T)

pscal=dt*sum(diag((phi_i)'*phi_j));