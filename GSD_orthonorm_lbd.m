function [lambda_coeffs_out] = GSD_orthonorm_lbd(PC, Kmodes, lambda_coeffs_in)

%Gram-Schmidt procedure:
ps11sq=pscal_lbd(PC, lambda_coeffs_in{1}(:), lambda_coeffs_in{1}(:));
lambda_coeffs_out{1}=lambda_coeffs_in{1}/sqrt(ps11sq);

for j=2:Kmodes
   tmp=lambda_coeffs_in{j};
   for k=1:j-1
       tmp=tmp-pscal_lbd(PC, lambda_coeffs_in{j}, lambda_coeffs_out{k})*lambda_coeffs_out{k};
   end
   %normalization:
   psjjsq=pscal_lbd(PC, tmp, tmp);
   lambda_coeffs_out{j}=tmp/sqrt(psjjsq);
end
