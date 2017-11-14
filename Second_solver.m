function [dis,vel,acc] = Second_solver(Method_TM2,u0,ut0,h,szt,K,C,M,F,show)
% %Newmark scheme
% %'Av', i.e. average scheme is unconditionally stable and Oh2 error
% %'FG', Fox&Goodwin has smallest error (Oh3) and conditionally stable (wh<2.45)
% tic;
reverseStr = '';
switch Method_TM2
    case 'Av'
        sigma =0.5;alpha = 0.25;
    case 'FG'
        sigma = 0.5;alpha = 1/12;
end

% %Newmark's scheme from Finite Element Procedures by Bathe
a0 = 1/(alpha*h^2);
a1 = sigma/(alpha*h);
a2 = 1/(alpha*h);
a3 = 1/(2*alpha)-1;
a4 = sigma/alpha-1;
a5 = h/2*(sigma/alpha-2);
a6 = h*(1-sigma);
a7 = sigma*h;
[LK,DK,PK] = ldl(K+a0*M+a1*C);%using LDL as recommend by Bathe
dis(:,1) = u0;
vel(:,1) = ut0;
utt0 = M\(F(:,1)-C*ut0-K*u0);
acc(:,1) = utt0;
if show
    fprintf('Time marching steps...');
end
opts.SYM = true;
for i = 2:szt
    load = F(:,i) + M*(a0*u0+a2*ut0+a3*utt0)+C*(a1*u0+a4*ut0+a5*utt0);
    u1 = PK*(LK'\(DK\(LK\(PK'*load))));
    dis(:,i) = u1;
    utt1 = a0*(u1-u0)-a2*ut0-a3*utt0;
    acc(:,i) = utt1;
    ut0 = ut0+a6*utt0+a7*utt1;
    vel(:,i) = ut0;
    u0 = u1;
    utt0 = utt1;
    if show
        %monitor the process
        str = sprintf('%.0f',i/szt*100);
        str = strcat(str,'%%\n');
        fprintf([reverseStr,str]);
        reverseStr = repmat(sprintf('\b'),1,length(str)-2);
    end
end
if show
    fprintf('\n');
end




