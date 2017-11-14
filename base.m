%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [base1] = base(psibasis, p, i, Ptot, Px1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -> Par's routine (unchanged).

% Construct the 1-var index structure psibase1
% of elements 0<= sum_s psibasis1{k}(s) <= totOrderPol
% Note psibase1 has all zeros except for variable i.
% this is a temporary variable used to identify one variable 
% parts of total basis.

  
  psibase1 = zeros(Px1,p);
  psibase1(:,i) = psibasis{1};

  base1 = zeros(Px1,1);
  ii = 1;
  jj = 1;
  %% test if row jj of the total basis is identical to row ii of psibase1, if so we have found a
  %% one variable element in the total base and we need to record this index
  while( ii<=Px1 )
    test = 0;
    for kk=1:p
      test = test + abs(psibasis{p}(jj,kk)-psibase1(ii,kk));
    end;
    if(test==0)
        base1(ii,1) = jj;
        ii = ii + 1;
    end; %if    
    jj = jj + 1;
  end    %while
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
