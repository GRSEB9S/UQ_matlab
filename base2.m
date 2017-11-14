%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [base2] = base2(psibasis, p, i, j, Ptot, Px2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -> Par's routine (unchanged).

% base2(): Mapping of two variable function to the tot basis
% see explanation for base()

  psibase2 = cell(Px2,1);
  tmp = zeros(1,p);
  for kk = 1:Px2
    psibase2{kk} = tmp;
    psibase2{kk}(i) = psibasis{2}(kk,1);
    psibase2{kk}(j) = psibasis{2}(kk,2);
  end;
  base2 = zeros(Px2,1);
  %
  %Match psibase2 for solution "i" with psibasis
  ii = 1;
  jj = 1;
  while( ii <= Px2 )
    test = 0;
    for kk=1:p
      test = test + abs(psibasis{p}(jj,kk) - psibase2{ii}(kk));
    end;
    if(test==0)
        base2(ii,1) = jj;
        ii = ii + 1;
    end; %if    
    jj = jj + 1;
  end    %while
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%