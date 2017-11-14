function [uijprime] = second_order(PC_2, M, K, C, F, nRandomVarsUsed, ...
    T_snap, dt_euler, nint, u0, ut0, Method_TM2, remain, remove, dfix)
reverseStr = '';
PA = nRandomVarsUsed*(nRandomVarsUsed-1)/2;
n = 0;
for i1=1:nRandomVarsUsed-1
    for i2=i1+1:nRandomVarsUsed
        
        newM = transformnew2(M,i1,i2);
        newK = transformnew2(K,i1,i2);
        newC = transformnew2(C,i1,i2);
        newF = transformnew2(F,i1,i2);
        [~,~,solGS, ~, ~] = GS(PC_2, newM, newK, newC, newF, ...
            nint, 2, T_snap, dt_euler, u0,ut0,Method_TM2,remain,remove,dfix,0);
        for j=1:PC_2.Size
            uijprime{i1}{i2}{j}=solGS((j-1)*nint+1:j*nint,:);
        end
        %monitor the process
        n = n+1;
        str = sprintf('%.0f',n/PA*100);
        str = strcat(str,'%%\n');
        fprintf([reverseStr,str]);
        reverseStr = repmat(sprintf('\b'),1,length(str)-2);
    end; 
end
end

function newX = transformnew2(X,i1,i2)
newX = cell(3,1);
newX{1} = X{1};
newX{2} = X{i1+1};
newX{3} = X{i2+1};
end

