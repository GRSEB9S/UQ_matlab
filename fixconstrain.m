function newX = fixconstrain(X,VorD,remain,remove,dfix)
%postprocess the result and apply the fixed points
szt = size(X,2);
nint = size(X,1);
switch VorD
    case {'Vel', 'Var'}
        newX = zeros(nint+length(dfix),szt);
        newX(remain,:) = X;
    case 'Dis'
        newX = zeros(nint+length(dfix),szt);
        newX(remain,:) = X;
        newX(remove,:) = kron(ones(1,szt),double(dfix'));
    case 'Comb'
        newnint = nint/2;
        vel = X(1:newnint,:);
        dis = X(newnint+1:nint,:);
        newvel = zeros(newnint+length(dfix),szt);
        newvel(remain,:) = vel;
        newdis(remain,:) = dis;
        newdis(remove,:) = kron(ones(1,szt),double(dfix'));
        newX = cat(1,newvel,newdis);
    case 'CombVar'
        newnint = nint/2;
        vel = X(1:newnint,:);
        dis = X(newnint+1:nint,:);
        newvel = zeros(newnint+length(dfix),szt);
        newvel(remain,:) = vel;
        newdis = zeros(newnint+length(dfix),szt);
        newdis(remain,:) = dis;
        newX = cat(1,newvel,newdis);
end
end