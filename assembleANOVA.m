function [Sol,uprime] = assembleANOVA(nRandomVarsUsed,orderPolSol,u0,uprime,ANOVA_ord,n,coeff)
%Assemble the ANOVA-Galerkin solution
Ptot=factorial(orderPolSol+nRandomVarsUsed)/(factorial(orderPolSol)*factorial(nRandomVarsUsed));
Sol=cell(Ptot,1);

psibasis=psibasisconstruct(nRandomVarsUsed,orderPolSol);
%Px1 & psibasis1
Px1=factorial(orderPolSol+1)/(factorial(orderPolSol)*factorial(1));
psibasis1=zeros(Px1,nRandomVarsUsed);

%%assemble the mean value, i.e. Sol{1}
for i = 1:2*n
    Sol{1}(i,:)=coeff{i}(1,:).*u0(i,:);
end

% Add one variable components:
for ii=1:nRandomVarsUsed
    uprime{1}{ii}{1}=uprime{1}{ii}{1}-u0;
    for i = 1:2*n
        Sol{1}(i,:)=Sol{1}(i,:)+coeff{i}(ii+1,:).*uprime{1}{ii}{1}(i,:);
    end
end
k = ii+1;
ind2 = zeros(nRandomVarsUsed,nRandomVarsUsed);
if ANOVA_ord >= 2 % Add two variables components:
    Px2=factorial(orderPolSol+2)/(factorial(orderPolSol)*factorial(2));
    psibasis2=zeros(Px2,nRandomVarsUsed,nRandomVarsUsed);
    for ii=1:nRandomVarsUsed-1
        for jj=(ii+1):nRandomVarsUsed
            uprime{2}{ii}{jj}{1}=uprime{2}{ii}{jj}{1}-uprime{1}{ii}{1}-uprime{1}{jj}{1}-u0;
            k = k+1;
            ind2(ii,jj)=k;
            for i = 1:2*n
                Sol{1}(i,:)=Sol{1}(i,:)+coeff{i}(k,:).*uprime{2}{ii}{jj}{1}(i,:);
            end
        end
    end
end
ind3 = zeros(nRandomVarsUsed,nRandomVarsUsed,nRandomVarsUsed);
if ANOVA_ord >= 3 % Add three variables components:
    Px3=factorial(orderPolSol+3)/(factorial(orderPolSol)*factorial(3));
    psibasis3 = cell(nRandomVarsUsed,nRandomVarsUsed,nRandomVarsUsed);
    for ii=1:nRandomVarsUsed-2
        for jj=(ii+1):nRandomVarsUsed-1
            for kk=(jj+1):nRandomVarsUsed
                uprime{3}{ii}{jj}{kk}{1}=uprime{3}{ii}{jj}{kk}{1}...
                    -uprime{2}{ii}{jj}{1}-uprime{2}{jj}{kk}{1}-uprime{2}{ii}{kk}{1}...
                    -uprime{1}{ii}{1}-uprime{1}{jj}{1}-uprime{1}{kk}{1}-u0;
                k = k+1;
                ind3(ii,jj,kk) = k;
                for i = 1:2*n
                    Sol{1}(i,:)=Sol{1}(i,:)+coeff{i}(k,:).*uprime{3}{ii}{jj}{kk}{1}(i,:);
                end
            end
        end
    end
end

%%assemble higher order components, i.e. Sol{i}
for ii=1:nRandomVarsUsed
    %book keeping for one component functions
    base1=base(psibasis,nRandomVarsUsed,ii,Ptot,Px1);
    psibasis1(1:Px1,ii)=base1(1:Px1,1);
end
if ANOVA_ord >= 2 %Px2 & psibasis2
    for ii=1:nRandomVarsUsed
        for jj=(ii+1):nRandomVarsUsed
            %book keeping for two component functions
            b2=base2(psibasis,nRandomVarsUsed,ii,jj,Ptot,Px2);
            psibasis2(1:Px2,ii,jj)=b2(1:Px2,1);
        end
    end
end
if ANOVA_ord == 3 %Px3 & psibasis3
    for ii=1:nRandomVarsUsed
        for jj=(ii+1):nRandomVarsUsed
            for kk=(jj+1):nRandomVarsUsed
                psibasis3{ii,jj,kk}=zeros(Px3,1);
                %book keeping for three component functions
                b3=base3(psibasis,nRandomVarsUsed,ii,jj,kk,Ptot,Px3);
                psibasis3{ii,jj,kk}(1:Px3)=b3(1:Px3,1);
            end
        end
    end
end

%assemble Sol{i}
for i = 2:Ptot
    Sol{i}=0*Sol{1};
end
for ii=1:nRandomVarsUsed
    for i=2:Px1
        ind=psibasis1(i,ii); % get correct global index for first component funcs.
        for j = 1:2*n
            Sol{ind}(j,:)=Sol{ind}(j,:)+coeff{j}(ii+1,:).*uprime{1}{ii}{i}(j,:);% Sol is total-global PC solution
        end
    end
end
if ANOVA_ord >= 2
    for ii=1:nRandomVarsUsed-1
        for jj=(ii+1):nRandomVarsUsed
            for i=2:Px2
                ind=psibasis2(i,ii,jj);
                for j= 1:2*n
                    Sol{ind}(j,:)=Sol{ind}(j,:)+coeff{j}(ind2(ii,jj),:).*uprime{2}{ii}{jj}{i}(j,:);
                end
            end
            for i=2:Px1
                indii=psibasis1(i,ii);
                indjj=psibasis1(i,jj);
                for j= 1:2*n
                    Sol{indii}(j,:)=Sol{indii}(j,:)-coeff{j}(ind2(ii,jj),:).*uprime{1}{ii}{i}(j,:);
                    Sol{indjj}(j,:)=Sol{indjj}(j,:)-coeff{j}(ind2(ii,jj),:).*uprime{1}{jj}{i}(j,:);
                end
            end
        end
    end
end
if ANOVA_ord == 3
    for ii=1:nRandomVarsUsed
        for jj=(ii+1):nRandomVarsUsed
            for kk=(jj+1):nRandomVarsUsed
                for i=2:Px3
                    ind=psibasis3{ii,jj,kk}(i);
                    for j = 1:2*n
                        Sol{ind}(j,:)=Sol{ind}(j,:)+coeff{j}(ind3(ii,jj,kk),:).*uprime{3}{ii}{jj}{kk}{i}(j,:);
                    end
                end
                for i=2:Px2
                    indij=psibasis2(i,ii,jj);
                    indik=psibasis2(i,ii,kk);
                    indjk=psibasis2(i,jj,kk);
                    for j =1:2*n
                        Sol{indij}(j,:)=Sol{indij}(j,:)-coeff{j}(ind3(ii,jj,kk),:).*uprime{2}{ii}{jj}{i}(j,:);
                        Sol{indik}(j,:)=Sol{indik}(j,:)-coeff{j}(ind3(ii,jj,kk),:).*uprime{2}{ii}{kk}{i}(j,:);
                        Sol{indjk}(j,:)=Sol{indjk}(j,:)-coeff{j}(ind3(ii,jj,kk),:).*uprime{2}{jj}{kk}{i}(j,:);
                    end
                end
                for i=2:Px1
                    indii=psibasis1(i,ii);
                    indjj=psibasis1(i,jj);
                    indkk=psibasis1(i,kk);
                    for j= 1:2*n
                        Sol{indii}(j,:)=Sol{indii}(j,:)+coeff{j}(ind3(ii,jj,kk),:).*uprime{1}{ii}{i}(j,:);
                        Sol{indjj}(j,:)=Sol{indjj}(j,:)+coeff{j}(ind3(ii,jj,kk),:).*uprime{1}{jj}{i}(j,:);
                        Sol{indkk}(j,:)=Sol{indkk}(j,:)+coeff{j}(ind3(ii,jj,kk),:).*uprime{1}{kk}{i}(j,:);
                    end
                end
            end
        end
    end
end

end