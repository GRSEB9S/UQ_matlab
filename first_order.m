function [uiprime] = first_order(PC_1, M, K, C, F, nRandomVarsUsed, ...
    T_snap, dt_euler, nint, u0, ut0, Method_TM2, remain, remove, dfix)
reverseStr = '';
n=0;
for i=1:nRandomVarsUsed
    newM = transformnew1(M,i);
    newK = transformnew1(K,i);
    newC = transformnew1(C,i);
    newF = transformnew1(F,i);
    [~,~,solGS,~, ~] = GS(PC_1, newM, newK, newC, newF, ...
        nint, 1, T_snap, dt_euler, u0, ut0,Method_TM2,remain,remove,dfix,0);
    for j=1:PC_1.Size
        uiprime{i}{j}=solGS((j-1)*nint+1:j*nint,:);
    end
    %monitor the process
    n = n+1;
    str = sprintf('%.0f',n/nRandomVarsUsed*100);
    str = strcat(str,'%%\n');
    fprintf([reverseStr,str]);
    reverseStr = repmat(sprintf('\b'),1,length(str)-2);
end
end

function newX = transformnew1(X,i)
newX = cell(2,1);
newX{1} = X{1};
newX{2} = X{i+1};
end
