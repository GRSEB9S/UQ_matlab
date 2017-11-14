%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [base3] = base3(psibasis, p, i, j, k, Ptot, Px3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -> Par's routine (unchanged).

% base2(): Mapping of two variable function to the tot basis
% see explanation for base()

  psibase3 = cell(Px3,1);
  tmp = zeros(1,p);
  for ll = 1:Px3
    psibase3{ll} = tmp;
    psibase3{ll}(i) = psibasis{3}(ll,1);
    psibase3{ll}(j) = psibasis{3}(ll,2);
    psibase3{ll}(k) = psibasis{3}(ll,3);
  end;
  base3 = zeros(Px3,1);
  %
  %Match psibase3 for solution "i" with psibasis
  ii = 1;
  jj = 1;
  while( ii <= Px3 )
    test = 0;
    for kk=1:p
      test = test + abs(psibasis{p}(jj,kk) - psibase3{ii}(kk));
    end;
    if(test==0)
        base3(ii,1) = jj;
        ii = ii + 1;
    end; %if    
    jj = jj + 1;
  end    %while
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%