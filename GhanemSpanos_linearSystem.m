%
% Implementation of the Ghanem Spanos PC Projection scheme
% to solve a stochastic linear system (Prasanth's code).
%

 function [xmean, xstd, xcell] = GhanemSpanos_linearSystem(A, b1, Psi2Norm, Dijk, P)

% A{i} 
% b1{}
% P - number of terms in the PC expansion for solution

  n = length(b1{1});

  Pa = size(A,1);
 
  BigA = sparse(zeros(n*P,n*P));
  
  BigA = kron(Dijk{1}(1:P,1:P),A{1});
  for i = 2:Pa
      temp = kron(Dijk{i}(1:P,1:P),A{i});
      BigA = BigA + temp;
  end 
  size(BigA);

  Pb = size(b1,1);

  Bigb = zeros(n*P,1);

  for i = 1:Pb
      Bigb(1+(i-1)*n:i*n,1) = Psi2Norm(i)*b1{i};
  end
  size(Bigb);

  x = BigA \ Bigb;

  x = reshape(x, n, P);

  xcell = cell(P,1);
  for i = 1:P
      xcell{i} = x(:,i);
  end 

  xmean=xcell{1}; %xmean = x(:,1); 

  covar = zeros(n,n);
  for i = 2:P
%     covar = covar + x(:,i)*x(:,i)'*Psi2Norm(i);
      covar = covar + xcell{i}*xcell{i}'*Psi2Norm(i);
  end
  xstd = sqrt(diag(covar));

% Compute true residual norm (not useful here)

%   Ax = PC_matvec(A, xcell, Psi2Norm, Dijk, P+1);
%  
%   r = cell(P+1,1);
%   for i = 1: Pb
%       r{i} = b1{i} - Ax{i};
%   end 
%   for i = Pb+1:P+1
%       r{i} = -Ax{i};
%   end
%  
%   res = 0;
%   for i = 1: P+1
%       res = res + r{i}'*r{i}*Psi2Norm(i);
%   end
%   res = sqrt(res);
