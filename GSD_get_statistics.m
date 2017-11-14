function [muGSD, varGSD] = GSD_get_statistics(phi_coeffs, lambda_coeffs, PC, Kmodes, nint, Nt)

%---compute the GSD mean:
muGSD=zeros(nint,Nt);
for i=1:Kmodes
    muGSD=muGSD+lambda_coeffs{i}(1)*phi_coeffs{i};
end

%---compute the GSD variance:
comb_lbdu=cell(PC.Size);
for p=1:PC.Size
   tmp=zeros(nint,Nt);
   for i=1:Kmodes
      tmp=tmp+lambda_coeffs{i}(p)*phi_coeffs{i};
   end
   comb_lbdu{p}=tmp;
end
varGSD=zeros(nint,Nt);
for p=1:PC.Size
   varGSD=varGSD+PC.PsiSqNorm(p)*(comb_lbdu{p}).^2;
end
varGSD=varGSD-PC.PsiSqNorm(1)*(comb_lbdu{1}).^2;