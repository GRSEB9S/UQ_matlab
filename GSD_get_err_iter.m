%
% compute error between 2 successive GSD expansions.
%

function [err_iter] = GSD_get_err_iter(phi_coeffs_new, lambda_coeffs_new, phi_coeffs_old, lambda_coeffs_old, ...
PC, Kmodes, nint, Nt, dt_euler)

%assemble coefficients of the difference:
for j=1:PC.Size
    tmp=zeros(nint,Nt);
    for i=1:Kmodes
       tmp=tmp+lambda_coeffs_new{i}(j)*phi_coeffs_new{i}-lambda_coeffs_old{i}(j)*phi_coeffs_old{i};
    end
    gamma{j}=tmp;
end

% for j=1:PC.Size
%     Mgamma{j}=massmat*gamma{j};
% end

err_iter=0.0;
for j=1:PC.Size
    err_iter=err_iter+PC.PsiSqNorm(j)*dt_euler*sum(diag((gamma{j})'*gamma{j}));
end



