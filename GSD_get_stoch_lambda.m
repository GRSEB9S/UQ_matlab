function PC_lambda=GSD_get_stoch_lambda(PC, lambda_coeffs, Kmodes, nRandomVarsUsed)

%size:
PC_lambda.Size=Kmodes;
%square of norms:
PC_lambda.PsiSqNorm=zeros(PC_lambda.Size,1);
matPC_sqnorm=diag(PC.PsiSqNorm);
for i=1:PC_lambda.Size
    PC_lambda.PsiSqNorm(i)=lambda_coeffs{i}'*matPC_sqnorm*lambda_coeffs{i};
end
%< \xi_m \lambda_i \lambda_j >:
PC_lambda.DD=cell(nRandomVarsUsed+1,1);
for m=1:nRandomVarsUsed+1
   PC_lambda.DD{m}=sparse(PC_lambda.Size,PC_lambda.Size);
   for i=1:PC_lambda.Size
     for j=1:PC_lambda.Size
        PC_lambda.DD{m}(i,j)=lambda_coeffs{i}'*PC.DD{m}*lambda_coeffs{j};
     end
   end
end
%< \xi_m \lambda_j >:
PC_lambda.Psilambda=zeros(nRandomVarsUsed+1,Kmodes);
for j=1:Kmodes
   for m=1:nRandomVarsUsed+1
        PC_lambda.Psilambda(m,j)=lambda_coeffs{j}(m)*PC.PsiSqNorm(m);
   end
end
%<\lambda_i>:
PC_lambda.MeanBasis=zeros(PC_lambda.Size,1);
for i=1:PC_lambda.Size
   PC_lambda.MeanBasis(i)=lambda_coeffs{i}(1);
end