function [mean_dis,std_dis] = MCS(M,K,C,F,u0,ut0,szt,dt,Method_TM2,totnRandomVarsUsed,xisample,remain,remove,dfix,N,output)
id = output.id;
reverseStr = '';
dis = zeros(N,szt);
for n = 1:N
    %assemble the coefficient matrices
    Msample = assemble(xisample(:,n),M,totnRandomVarsUsed);
    Ksample = assemble(xisample(:,n),K,totnRandomVarsUsed);
    Csample = assemble(xisample(:,n),C,totnRandomVarsUsed);
    Fsample = assemble(xisample(:,n),F,totnRandomVarsUsed);
    
    %second order solver
    [stockdis,~,~] = Second_solver(Method_TM2,u0,ut0,dt,szt,Ksample,Csample,Msample,Fsample,0);
    stockdis = fixconstrain(stockdis,'Dis',remain,remove,dfix);
    
    %output only the indicated nodes
    dis(n,:) = stockdis(id,:);
    
    %monitor the process
    str = sprintf('%.0f',n/N*100);
    str = strcat(str,'%%\n');
    fprintf([reverseStr,str]);
    reverseStr = repmat(sprintf('\b'),1,length(str)-2);
end
mean_dis = mean(dis);
std_dis = std(dis);

function Xassembled = assemble(xisample,X,n)
Xassembled = X{1};
for i = 1:n
    Xassembled = Xassembled + X{i+1}*xisample(i);
end 