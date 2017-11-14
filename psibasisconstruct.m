%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [psibasis] =  psibasisconstruct(nRandVar,order)   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ->Par's routine (unchanged).
% Modified by Lin Gao , lingao.gao@mail.utoronto.ca
% Aug 17, 2015

psibasis=cell(nRandVar,1);
for i=1:nRandVar
    basis=BASIS(i,order);
    P = int16(factorial(i+order)/(factorial(i)*factorial(order)));
    for j=1:P
      for k=1:i
        psibasis{i}(j,k)=basis(j,i+1-k);
      end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
